// SPDX-License-Identifier: MIT

/**
 * File: gpu_discretizer.h
 * Project: foldcomp
 * Created: 2025-11-05
 * Description:
 *     GPU-accelerated discretization operations
 */

#pragma once
#include "gpu_buffer.h"
#include <cstddef>
#include <cstdint>
#include <cuda_runtime.h>
#include <memory>
#include <vector>

/**
 * @brief GPU-accelerated context for discretization operations
 *
 * This class manages GPU resources and provides methods for batch
 * continuization of discretized values using CUDA.
 */
class GPUDiscretizerContext {
public:
    GPUDiscretizerContext();
    ~GPUDiscretizerContext();

    // Disable copy and move
    GPUDiscretizerContext(const GPUDiscretizerContext &) = delete;
    GPUDiscretizerContext &operator=(const GPUDiscretizerContext &) = delete;
    GPUDiscretizerContext(GPUDiscretizerContext &&) = delete;
    GPUDiscretizerContext &operator=(GPUDiscretizerContext &&) = delete;

    /**
   * @brief Continuize backbone angles, temperature factors, and sidechain
   * angles in a single GPU round trip (one sync instead of three).
   *
   * Uploads all three input arrays, launches all three kernels, enqueues any
   * requested D→H transfers — but does NOT sync. All enqueued work is
   * stream-ordered before the caller's final cudaEventSynchronize on its
   * per-slot CPU-sync event (see fillAtomCoordinatesGPU), so a mid-pipeline
   * sync here is unnecessary.
   */
    void continuize_all(
        const uint64_t *bb_raw, // [n_bb_residues] BackboneChain as uint64
        int n_bb_residues,
        const std::vector<float>
            &bb_mins_all, // [n_structs * 6] per-struct min values
        const std::vector<float>
            &bb_cont_fs_all, // [n_structs * 6] per-struct cont_f values
        const std::vector<int>
            &bb_struct_offsets,  // [n_structs+1] cumulative residue counts
        PinnedFloatVec &bb_cont, // pinned; unused if out_d_bb != nullptr
        const std::vector<uint8_t> &tf_disc, const std::vector<int> &tf_offsets,
        const std::vector<float> &tf_mins, const std::vector<float> &tf_cont_fs,
        const std::vector<uint8_t, CudaPinnedAllocator<uint8_t>> &sc_disc,
        float sc_min, float sc_cont_f,
        cudaStream_t stream,        // compute stream for this pipeline
        float **out_d_bb = nullptr, // backbone stays on device
        float **out_h_tf = nullptr, // TF pinned ptr (no heap copy)
        float **out_d_sc = nullptr, // sidechain stays on device (no D→H)
        float **out_d_tf =
            nullptr // TF device ptr (skip H→D re-upload in fill kernel)
    );

    /**
   * @brief Enable or disable fused kernels (default: enabled)
   *
   * When enabled, each continuize operation uses a single fused kernel that
   * reads uint/uint8_t and writes float in one pass, halving global memory
   * traffic.
   */
    void setUseFused(bool v);

    /**
   * @brief Pre-allocate all internal device/pinned buffers to worst-case sizes.
   * Call once at startup (before the first batch) to avoid cold-start
   * cudaMalloc latency.
   *
   * @param max_batch_size  Maximum number of structures per batch
   * @param max_residues    Upper bound on residues per structure
   * @param max_sc_per_res  Upper bound on sidechain torsion angles per residue
   * (default 4)
   */
    void preallocate(int max_batch_size, int max_residues, int max_sc_per_res = 4,
        cudaStream_t stream = 0);

private:
    struct Impl;
    std::unique_ptr <Impl> impl;
};
