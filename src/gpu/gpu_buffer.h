// SPDX-License-Identifier: MIT
/**
 * File: gpu_buffer.h
 * Project: foldcomp
 * Description:
 *     Grow-only GPU device buffer helper.
 *     Use as a static or long-lived object to eliminate per-call cudaMalloc overhead.
 *     The buffer grows when a larger allocation is requested but never shrinks.
 *     Suitable for persistent per-function or per-module device memory.
 */
#pragma once
#include <cuda_runtime.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>

#include "gpu_stream.h"

/**
 * Grow-only pinned (page-locked) host buffer.
 * Use for staging buffers that participate in cudaMemcpyAsync to/from device.
 * Pinned memory enables true DMA without CUDA's internal bounce-buffer copy.
 */
struct GrowOnlyPinnedBuf {
    void* ptr = nullptr;
    size_t capacity = 0;

    void* ensure(size_t bytes) {
        if (bytes > capacity) {
            size_t new_cap = std::max(bytes, capacity * 2);
            CUDA_CHECK_FREE(cudaFreeHost(ptr));
            ptr = nullptr;
            capacity = 0;
            if (cudaMallocHost(&ptr, new_cap) != cudaSuccess)
                throw std::runtime_error("GrowOnlyPinnedBuf: cudaMallocHost failed");
            capacity = new_cap;
        }
        return ptr;
    }

    template <typename T>
    T* ensure_n(size_t n) { return static_cast<T*>(ensure(n * sizeof(T))); }

    /** Pre-allocate at least `bytes` bytes. Call at startup to avoid cold-start cudaMallocHost latency. */
    void preallocate(size_t bytes) { ensure(bytes); }

    ~GrowOnlyPinnedBuf() { cudaFreeHostChecked(ptr); }

    GrowOnlyPinnedBuf() = default;
    GrowOnlyPinnedBuf(const GrowOnlyPinnedBuf&) = delete;
    GrowOnlyPinnedBuf& operator=(const GrowOnlyPinnedBuf&) = delete;
};

struct GrowOnlyDeviceBuf {
    void* ptr = nullptr;
    size_t capacity = 0;

    /**
     * Ensure at least `bytes` bytes are allocated on the device.
     * Returns the device pointer. Reallocates only when capacity is insufficient.
     *
     * Uses cudaMallocFromPoolAsync/cudaFreeAsync (stream-ordered, pool-backed) against
     * Foldcomp's private per-device pool (getFoldcompMemPool, gpu_stream.h) rather than
     * cudaMalloc/cudaFree or the device's default pool. Freed memory is cached and
     * reused on regrow within that private pool, without holding it in this buffer's
     * own free-list and without making the pool's release-threshold tuning (kept
     * permanently high to avoid reallocation latency) visible to or shared with any
     * other library's allocations on the same device.
     */
    void* ensure(size_t bytes, cudaStream_t stream = 0) {
        if (bytes > capacity) {
            size_t new_cap = std::max(bytes, capacity * 2);
            CUDA_CHECK_FREE(cudaFreeAsync(ptr, stream));
            ptr = nullptr;
            capacity = 0;
            if (cudaMallocFromPoolAsync(&ptr, new_cap, getFoldcompMemPool(), stream) != cudaSuccess)
                throw std::runtime_error("GrowOnlyDeviceBuf: cudaMallocFromPoolAsync failed");
            capacity = new_cap;
        }
        return ptr;
    }

    /** Typed convenience wrapper: ensure n elements of type T. */
    template <typename T>
    T* ensure_n(size_t n, cudaStream_t stream = 0) { return static_cast<T*>(ensure(n * sizeof(T), stream)); }

    /** Pre-allocate at least `bytes` bytes. Call at startup to avoid cold-start cudaMallocAsync latency. */
    void preallocate(size_t bytes, cudaStream_t stream = 0) { ensure(bytes, stream); }

    // Plain synchronous cudaFree (not cudaFreeAsync) is intentional here: it is a
    // documented valid way to release memory obtained from cudaMallocAsync, it
    // still returns the memory to the shared pool, and it needs no stream handle
    // — the owning stream may already be destroyed by the time this destructor
    // runs (e.g. FoldcompGPUState tears down its streams before its buffer
    // members destruct).
    ~GrowOnlyDeviceBuf() { cudaFreeChecked(ptr); }

    GrowOnlyDeviceBuf() = default;
    GrowOnlyDeviceBuf(const GrowOnlyDeviceBuf&) = delete;
    GrowOnlyDeviceBuf& operator=(const GrowOnlyDeviceBuf&) = delete;
};

/**
 * STL-compatible allocator backed by CUDA pinned (page-locked) memory.
 * Lets std::vector serve as a direct DMA target for cudaMemcpyAsync,
 * eliminating the intermediate staging memcpy.
 */
template <typename T>
struct CudaPinnedAllocator {
    using value_type = T;
    T* allocate(std::size_t n) {
        T* ptr = nullptr;
        if (cudaMallocHost(&ptr, n * sizeof(T)) != cudaSuccess)
            throw std::bad_alloc();
        return ptr;
    }
    void deallocate(T* ptr, std::size_t) noexcept { cudaFreeHostChecked(ptr); }
    template <typename U>
    bool operator==(const CudaPinnedAllocator<U>&) const noexcept { return true; }
    template <typename U>
    bool operator!=(const CudaPinnedAllocator<U>&) const noexcept { return false; }
};

using PinnedFloatVec = std::vector<float, CudaPinnedAllocator<float>>;

/** Round `off` up to the next multiple of `align` (must be power of 2). */
static inline size_t gpu_align_up(size_t off, size_t align) {
    return (off + align - 1) & ~(align - 1);
}
