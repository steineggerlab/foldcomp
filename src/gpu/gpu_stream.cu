// SPDX-License-Identifier: MIT

#include "gpu_stream.h"
#include <mutex>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>

// Internal streams used only to create/pre-record the CPU-sync event pool below
// (see initGPUPipelineEvents/destroyGPUPipelineEvents). Not exposed outside this
// file — every FoldcompGPUState owns its own GPUStreamEventPipeline instead.
// CUDA streams/events are per-device, so one pipeline is kept per device, keyed
// by device id (same pattern as gpu_decompression_pipeline.cpp's device_id-keyed
// pipeline pool).
static std::unordered_map<int, GPUStreamEventPipeline> g_pipes;

// Pool of per-batch CPU-sync events, shared across every pooled GPU pipeline instance
// on the same device (see gpu_decompression_pipeline.cpp's PooledPipeline pool).
// Events are fungible within a device — any event works identically for any caller
// on that device — but not across devices, so the pool is keyed by device id, one
// grow-only queue per device sized to the sum of gpu_slots across live pool
// instances on that device. Guarded by its own mutex so concurrently-running
// pipeline instances can safely acquire/release without relying on a single global
// run-lock.
static std::mutex g_cpu_sync_mutex;
static std::unordered_map<int, std::queue<cudaEvent_t>> g_cpu_sync_pools;

// Foldcomp's private per-device cudaMallocAsync pool (see getFoldcompMemPool in
// gpu_stream.h). Kept separate from g_cpu_sync_mutex since pool lookups happen on
// every GrowOnlyDeviceBuf growth, independent of the event-pool machinery above.
static std::mutex g_mem_pool_mutex;
static std::unordered_map<int, cudaMemPool_t> g_mem_pools;

void createGPUStreamEventPipeline(GPUStreamEventPipeline& pipe) {
    CUDA_CHECK(cudaStreamCreate(&pipe.stream_a));
    CUDA_CHECK(cudaStreamCreate(&pipe.stream_b));
    CUDA_CHECK(cudaEventCreateWithFlags(&pipe.fill_done, cudaEventDisableTiming));
    CUDA_CHECK(cudaEventCreateWithFlags(&pipe.slot_ready, cudaEventDisableTiming));
    // Pre-record as "completed" so the first fill kernel needs no guard.
    CUDA_CHECK(cudaEventRecord(pipe.slot_ready, pipe.stream_a));
}

void destroyGPUStreamEventPipeline(GPUStreamEventPipeline& pipe) {
    if (pipe.stream_a)  { cudaStreamSynchronize(pipe.stream_a);  cudaStreamDestroy(pipe.stream_a);  pipe.stream_a  = nullptr; }
    if (pipe.stream_b)  { cudaStreamSynchronize(pipe.stream_b);  cudaStreamDestroy(pipe.stream_b);  pipe.stream_b  = nullptr; }
    if (pipe.fill_done)  { cudaEventDestroy(pipe.fill_done);  pipe.fill_done  = nullptr; }
    if (pipe.slot_ready) { cudaEventDestroy(pipe.slot_ready); pipe.slot_ready = nullptr; }
}

cudaEvent_t acquireCPUSyncEvent() {
    int device = 0;
    cudaGetDevice(&device);
    std::lock_guard<std::mutex> lock(g_cpu_sync_mutex);
    std::queue<cudaEvent_t>& pool = g_cpu_sync_pools[device];
    if (pool.empty())
        throw std::runtime_error("CPU sync event pool exhausted");
    cudaEvent_t e = pool.front();
    pool.pop();
    return e;
}

void releaseCPUSyncEvent(cudaEvent_t e) {
    int device = 0;
    cudaGetDevice(&device);
    std::lock_guard<std::mutex> lock(g_cpu_sync_mutex);
    g_cpu_sync_pools[device].push(e);
}

void initGPUPipelineEvents(int n_in_flight) {
    int device = 0;
    cudaGetDevice(&device);
    std::lock_guard<std::mutex> lock(g_cpu_sync_mutex);
    GPUStreamEventPipeline& pipe = g_pipes[device];
    if (!pipe.stream_a && !pipe.stream_b && !pipe.fill_done && !pipe.slot_ready) {
        createGPUStreamEventPipeline(pipe);
    }
    // Per-batch CPU sync events: one per in-flight slot, aggregated across every
    // pipeline instance resident on this device (see comment above g_cpu_sync_pools).
    // Pre-record on this device's stream_a so they start in "completed" state.
    std::queue<cudaEvent_t>& pool = g_cpu_sync_pools[device];
    while ((int)pool.size() < n_in_flight) {
        cudaEvent_t cpu_evt = nullptr;
        CUDA_CHECK(cudaEventCreateWithFlags(&cpu_evt, cudaEventDisableTiming));
        CUDA_CHECK(cudaEventRecord(cpu_evt, pipe.stream_a));
        pool.push(cpu_evt);
    }
}

cudaMemPool_t getFoldcompMemPool() {
    int device = 0;
    cudaGetDevice(&device);
    std::lock_guard<std::mutex> lock(g_mem_pool_mutex);
    cudaMemPool_t& pool = g_mem_pools[device];
    if (!pool) {
        cudaMemPoolProps props{};
        props.allocType = cudaMemAllocationTypePinned;
        props.location.type = cudaMemLocationTypeDevice;
        props.location.id = device;
        CUDA_CHECK(cudaMemPoolCreate(&pool, &props));
        // Keep freed blocks cached in this pool across batches (avoids reintroducing
        // allocation latency on the next batch-size ramp-up). Safe to pin at max here
        // — unlike the device's default pool, this pool is private to Foldcomp, so it
        // has no effect on any other library's cudaMallocAsync calls on this device.
        uint64_t threshold = UINT64_MAX;
        CUDA_CHECK(cudaMemPoolSetAttribute(pool, cudaMemPoolAttrReleaseThreshold, &threshold));
    }
    return pool;
}

void destroyGPUPipelineEvents() {
    std::lock_guard<std::mutex> lock(g_cpu_sync_mutex);
    int active_device = 0;
    cudaGetDevice(&active_device);
    for (auto& kv : g_pipes) {
        cudaSetDevice(kv.first);
        destroyGPUStreamEventPipeline(kv.second);
    }
    g_pipes.clear();
    for (auto& kv : g_cpu_sync_pools) {
        cudaSetDevice(kv.first);
        while (!kv.second.empty()) {
            cudaEventDestroy(kv.second.front());
            kv.second.pop();
        }
    }
    g_cpu_sync_pools.clear();
    {
        std::lock_guard<std::mutex> pool_lock(g_mem_pool_mutex);
        for (auto& kv : g_mem_pools) {
            cudaSetDevice(kv.first);
            if (kv.second) cudaMemPoolDestroy(kv.second);
        }
        g_mem_pools.clear();
    }
    cudaSetDevice(active_device);
}
