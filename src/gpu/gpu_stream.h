// SPDX-License-Identifier: MIT

#pragma once

#include <cstdio>
#include <cuda_runtime.h>

// Throw a std::runtime_error on any non-success CUDA return value.
// Callers must include <stdexcept> and <string> (or include gpu_stream.h from a .cu/.cpp that does).
// Usage: CUDA_CHECK(cudaMemcpyAsync(...));
#define CUDA_CHECK(call) \
    do { \
        cudaError_t _e = (call); \
        if (_e != cudaSuccess) \
            throw std::runtime_error( \
                std::string(__FILE__) + ":" + std::to_string(__LINE__) + \
                " CUDA error: " + cudaGetErrorString(_e)); \
    } while (0)

// Like CUDA_CHECK, but for cudaFree/cudaFreeHost calls: tolerates
// cudaErrorCudartUnloading, which the driver returns instead of doing the
// free when the CUDA runtime has already begun tearing down (e.g. a
// free reached from a static/global destructor after main() returns).
// That is an expected no-op, not a real failure — see DALI's handling of
// the same error for the same reason. Any other non-success value still
// throws. Callers must include <stdexcept> and <string>.
#define CUDA_CHECK_FREE(call) \
    do { \
        cudaError_t _e = (call); \
        if (_e != cudaSuccess && _e != cudaErrorCudartUnloading) \
            throw std::runtime_error( \
                std::string(__FILE__) + ":" + std::to_string(__LINE__) + \
                " CUDA error: " + cudaGetErrorString(_e)); \
    } while (0)

// noexcept variants for use in destructors, where throwing is not an option.
// Unexpected errors are reported to stderr instead of being thrown;
// cudaErrorCudartUnloading is silently ignored (see CUDA_CHECK_FREE above).
inline void cudaFreeChecked(void* ptr) noexcept {
    cudaError_t e = cudaFree(ptr);
    if (e != cudaSuccess && e != cudaErrorCudartUnloading)
        std::fprintf(stderr, "cudaFree failed: %s\n", cudaGetErrorString(e));
}

inline void cudaFreeHostChecked(void* ptr) noexcept {
    cudaError_t e = cudaFreeHost(ptr);
    if (e != cudaSuccess && e != cudaErrorCudartUnloading)
        std::fprintf(stderr, "cudaFreeHost failed: %s\n", cudaGetErrorString(e));
}

/**
 * Number of independent GPU slots.  Each slot owns its own stream pair (stream_a/stream_b),
 * single output buffer, and persistent GPU buffers.  MAX_IN_FLIGHT batches can be in-flight
 * simultaneously: while slot K is doing D→H on stream_b, slots K+1…K+N-1 run kernels on
 * their stream_a, hiding D→H latency without any intra-slot rotation bookkeeping.
 */
static constexpr int MAX_IN_FLIGHT = 4;

/**
 * A two-stream pipeline: stream_a runs kernels, stream_b handles D→H transfers,
 * and fill_done/slot_ready synchronize the two without a device-wide sync.
 * Every FoldcompGPUState owns one of these (gpu_foldcomp_batch.cpp), created via
 * initOwnedStreams()/destroyOwnedStreams(), which just forward to the two
 * functions below so the create/destroy logic isn't duplicated per slot.
 */
struct GPUStreamEventPipeline {
    cudaStream_t stream_a  = nullptr;
    cudaStream_t stream_b  = nullptr;
    cudaEvent_t  fill_done  = nullptr;
    cudaEvent_t  slot_ready = nullptr;
};

/**
 * Creates stream_a/stream_b and fill_done/slot_ready events (all fields must
 * be nullptr on entry), and pre-records slot_ready on stream_a as "completed"
 * so the first fill kernel needs no guard.
 */
void createGPUStreamEventPipeline(GPUStreamEventPipeline& pipe);

/**
 * Synchronizes and destroys every non-null field of pipe, resetting them to
 * nullptr. Safe to call on a partially- or fully-nulled pipe.
 */
void destroyGPUStreamEventPipeline(GPUStreamEventPipeline& pipe);

/**
 * Per-batch CPU-sync event pool, keyed internally by the current device (acquire/init
 * operate on whichever device is current via cudaGetDevice, mirroring the device_id-keyed
 * pipeline pool in gpu_decompression_pipeline.cpp). Each in-flight batch holds one event
 * that is recorded on stream B after D→H; the CPU calls cudaEventSynchronize on it before
 * consuming results. Events are only ever returned to the pool after synchronization, so
 * re-recording is safe. An event acquired on one device must be released while that same
 * device is current — events are not fungible across devices.
 *
 * Pool size == MAX_IN_FLIGHT (one per in-flight slot), per device.
 */
cudaEvent_t acquireCPUSyncEvent();
void releaseCPUSyncEvent(cudaEvent_t e);

/**
 * Grows the current device's cpu-sync event pool to hold at least n_in_flight events
 * (pre-recorded as "already done" so the first user of each needs no guard). Grow-only:
 * never shrinks. Safe to call repeatedly as pipeline instances are built.
 * @param n_in_flight  Minimum pool size for the current device (default: MAX_IN_FLIGHT).
 */
void initGPUPipelineEvents(int n_in_flight = MAX_IN_FLIGHT);

/** Destroy all pipeline events, drain the cpu-sync pools, and destroy the private
 * memory pools (see getFoldcompMemPool below) for every device.  Call at program exit. */
void destroyGPUPipelineEvents();

/**
 * Foldcomp's private, per-device cudaMallocAsync memory pool (get-or-create, keyed by
 * the current device via cudaGetDevice — same pattern as acquireCPUSyncEvent above).
 * GrowOnlyDeviceBuf allocates from this pool instead of the device's default pool so
 * that raising its release threshold (to keep freed blocks cached across batches
 * instead of eagerly trimming them back to the OS) only affects Foldcomp's own
 * allocations, never a downstream library's cudaMallocAsync calls sharing the device.
 */
cudaMemPool_t getFoldcompMemPool();
