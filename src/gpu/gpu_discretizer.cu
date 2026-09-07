// SPDX-License-Identifier: MIT

/**
 * File: gpu_discretizer.cu
 * Project: foldcomp
 * Created: 2025-11-05
 * Description:
 *     GPU-accelerated discretization operations using
 */

#include "gpu_discretizer.h"
#include "gpu_buffer.h"
#include "gpu_stream.h"
#include <cuda_runtime.h>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <nvtx3/nvtx3.hpp>

// Per-struct backbone kernel: unpacks raw FCZ bytes (8 bytes/residue) + continuizes in one pass.
// Binary-searches bb_offsets[0..n_structs] to find which struct residue i belongs to,
// then applies that struct's mins/cont_fs (6 values each).
__global__ void unpack_raw_file_continuize_bb_per_struct_kernel(
    const uint64_t* __restrict__ bb_raw,  // [n_residues]
    float* __restrict__ out,              // [n_residues * 6]
    const int* __restrict__ bb_offsets,   // [n_structs+1] cumulative residue counts
    const float* __restrict__ bb_mins,    // [n_structs * 6]
    const float* __restrict__ bb_cont_fs, // [n_structs * 6]
    int n_structs,
    int n_residues) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_residues) return;

    // Binary search for structure index
    int lo = 0, hi = n_structs - 1;
    while (lo < hi) {
        int mid = (lo + hi) >> 1;
        if (bb_offsets[mid + 1] <= i)
            lo = mid + 1;
        else
            hi = mid;
    }
    const float* mins = bb_mins + lo * 6;
    const float* cont_fs = bb_cont_fs + lo * 6;

    uint64_t v = bb_raw[i];
    uint32_t omega = (uint32_t)((v & 0x07u) << 8) | (uint32_t)((v >> 8) & 0xFFu);
    uint32_t psi = (uint32_t)((v >> 12) & 0xFF0u) | (uint32_t)((v >> 28) & 0x00Fu);
    uint32_t phi = (uint32_t)((v >> 16) & 0x0F00u) | (uint32_t)((v >> 32) & 0xFFu);
    uint32_t ca_c_n = (uint32_t)((v >> 40) & 0xFFu);
    uint32_t c_n_ca = (uint32_t)((v >> 48) & 0xFFu);
    uint32_t n_ca_c = (uint32_t)((v >> 56) & 0xFFu);
    out[i * 6 + 0] = (float)phi * cont_fs[0] + mins[0];    // phi
    out[i * 6 + 1] = (float)psi * cont_fs[1] + mins[1];    // psi
    out[i * 6 + 2] = (float)omega * cont_fs[2] + mins[2];  // omega
    out[i * 6 + 3] = (float)n_ca_c * cont_fs[3] + mins[3]; // n_ca_c
    out[i * 6 + 4] = (float)ca_c_n * cont_fs[4] + mins[4]; // ca_c_n
    out[i * 6 + 5] = (float)c_n_ca * cont_fs[5] + mins[5]; // c_n_ca
}

// uint8_t variant of continuize_fused_kernel (SC angles: values 0-255)
__global__ void continuize_fused_u8_kernel(
    const uint8_t* input, float* output, int n, float min_v, float cont_f) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) output[i] = (float)input[i] * cont_f + min_v;
}

// uint8_t variant of continuize_per_struct_kernel (TF: per-struct min/cont_f)
__global__ void continuize_per_struct_u8_kernel(
    const uint8_t* d_input,
    float* d_output,
    const int* d_offsets, // n_structs+1 values
    const float* d_mins,
    const float* d_cont_fs,
    int n_structs,
    int total) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;
    int lo = 0, hi = n_structs - 1;
    while (lo < hi) {
        int mid = (lo + hi) >> 1;
        if (d_offsets[mid + 1] <= idx)
            lo = mid + 1;
        else
            hi = mid;
    }
    d_output[idx] = d_mins[lo] + (float)d_input[idx] * d_cont_fs[lo];
}

struct GPUDiscretizerContext::Impl {
    bool useFused;
    Impl() : useFused(true) {}
    ~Impl() {}

    // continuize_all buffers
    GrowOnlyDeviceBuf all_d_bb_raw, all_d_bb_out;
    GrowOnlyDeviceBuf all_d_tf_pack, all_d_tf_out;
    GrowOnlyDeviceBuf all_d_sc, all_d_sc_out;
    GrowOnlyDeviceBuf all_d_bb_params;
    GrowOnlyPinnedBuf all_ph_tf_pack, all_ph_tf_out;
    GrowOnlyPinnedBuf all_ph_sc_out;
    GrowOnlyPinnedBuf all_ph_bb_params;
};

GPUDiscretizerContext::GPUDiscretizerContext() : impl(new Impl()) {}

GPUDiscretizerContext::~GPUDiscretizerContext() { }

void GPUDiscretizerContext::setUseFused(bool v) { impl->useFused = v; }

void GPUDiscretizerContext::preallocate(int max_batch_size, int max_residues, int max_sc_per_res,
    cudaStream_t stream) {
    const int n_bb = max_batch_size * max_residues;
    const int n_tf = n_bb;
    const int n_sc = n_bb * max_sc_per_res;
    const int n_str = max_batch_size;

    impl->all_d_bb_raw.ensure_n<uint64_t>(n_bb, stream);
    impl->all_d_bb_out.ensure_n<float>(n_bb * 6, stream);
    impl->all_d_tf_out.ensure_n<float>(n_tf, stream);
    impl->all_d_sc.ensure_n<uint8_t>(n_sc, stream);
    impl->all_d_sc_out.ensure_n<float>(n_sc, stream);

    // bb_params pack: (n_str+1)*int + n_str*6*float (mins) + n_str*6*float (cont_fs)
    const size_t bb_params_sz = (size_t)(n_str + 1) * sizeof(int) + (size_t)n_str * 12 * sizeof(float);
    impl->all_d_bb_params.ensure(bb_params_sz, stream);
    impl->all_ph_bb_params.ensure(bb_params_sz);

    // tf_pack: n_tf*u8 + align pad + (n_str+1)*int + n_str*float (mins) + n_str*float (cont_fs)
    const size_t tf_pack_sz = gpu_align_up((size_t)n_tf, 4) + (size_t)(n_str + 1) * sizeof(int) + (size_t)n_str * 2 * sizeof(float);
    impl->all_d_tf_pack.ensure(tf_pack_sz, stream);
    impl->all_ph_tf_pack.ensure(tf_pack_sz);
    impl->all_ph_tf_out.ensure_n<float>(n_tf);
}

/// Combined: upload all three arrays, launch all three kernels, enqueue D→H
// transfers — no internal sync. Caller syncs once in fillAtomCoordinatesGPU.
void GPUDiscretizerContext::continuize_all(
    const uint64_t* bb_raw,
    int n_bb_residues,
    const std::vector<float>& bb_mins_all,     // n_structs * 6 per-struct min values
    const std::vector<float>& bb_cont_fs_all,  // n_structs * 6 per-struct cont_f values
    const std::vector<int>& bb_struct_offsets, // n_structs+1 cumulative residue counts
    PinnedFloatVec& bb_cont,
    const std::vector<uint8_t>& tf_disc,
    const std::vector<int>& tf_offsets,
    const std::vector<float>& tf_mins,
    const std::vector<float>& tf_cont_fs,
    const std::vector<uint8_t, CudaPinnedAllocator<uint8_t>>& sc_disc,
    float sc_min,
    float sc_cont_f,
    cudaStream_t stream,
    float** out_d_bb,
    float** out_h_tf,
    float** out_d_sc,
    float** out_d_tf) {
    int n_bb_out = n_bb_residues * 6; // floats in continuized backbone output
    int n_tf = (int)tf_disc.size();
    int n_sc = (int)sc_disc.size();
    int n_structs = (int)tf_mins.size();

    // --- Device buffers ---
    auto* d_bb_raw = impl->all_d_bb_raw.ensure_n<uint64_t>(n_bb_residues, stream);
    auto* d_bb_out = impl->all_d_bb_out.ensure_n<float>(n_bb_out, stream);

    // TF pack offsets: [uint8_t[n_tf] tf | pad | int[n+1] off | float[n] min | float[n] cf]
    const size_t sz_tf_data = (size_t)n_tf * sizeof(uint8_t);
    const size_t off_tf_data = 0;
    const size_t off_tf_off = (n_tf > 0 && n_structs > 0) ? gpu_align_up(sz_tf_data, 4) : 0;
    const size_t sz_tf_off = (size_t)(n_structs + 1) * sizeof(int);
    const size_t off_tf_min = off_tf_off + sz_tf_off;
    const size_t sz_tf_min = (size_t)n_structs * sizeof(float);
    const size_t off_tf_cf = off_tf_min + sz_tf_min;
    const size_t sz_tf_cf = (size_t)n_structs * sizeof(float);
    const size_t tf_pack_sz = (n_tf > 0 && n_structs > 0) ? (off_tf_cf + sz_tf_cf) : 0;

    uint8_t* d_tf = nullptr;
    float* d_tf_out = nullptr;
    int* d_tf_off = nullptr;
    float* d_tf_min = nullptr;
    float* d_tf_cf = nullptr;
    if (tf_pack_sz > 0) {
        char* dp = static_cast<char*>(impl->all_d_tf_pack.ensure(tf_pack_sz, stream));
        d_tf = reinterpret_cast<uint8_t*>(dp + off_tf_data);
        d_tf_off = reinterpret_cast<int*>(dp + off_tf_off);
        d_tf_min = reinterpret_cast<float*>(dp + off_tf_min);
        d_tf_cf = reinterpret_cast<float*>(dp + off_tf_cf);
        d_tf_out = impl->all_d_tf_out.ensure_n<float>(n_tf, stream);
    }

    uint8_t* d_sc = nullptr;
    float* d_sc_out = nullptr;
    if (n_sc > 0) {
        d_sc = impl->all_d_sc.ensure_n<uint8_t>(n_sc, stream);
        d_sc_out = impl->all_d_sc_out.ensure_n<float>(n_sc, stream);
    }

    // --- Pinned staging (TF only; bb and sc are already pinned at call site) ---

    // Backbone: bb_raw is already pinned (compressedBackBone_flat uses CudaPinnedAllocator)
    CUDA_CHECK(cudaMemcpyAsync(d_bb_raw, bb_raw, n_bb_residues * sizeof(uint64_t), cudaMemcpyHostToDevice, stream));

    // TF: pack [data|off|min|cf] into a single pinned+device buffer, one H→D transfer
    float* h_tf_out = nullptr;
    if (tf_pack_sz > 0) {
        char* hp = static_cast<char*>(impl->all_ph_tf_pack.ensure(tf_pack_sz));
        h_tf_out = impl->all_ph_tf_out.ensure_n<float>(n_tf);
        memcpy(hp + off_tf_data, tf_disc.data(), sz_tf_data);
        memcpy(hp + off_tf_off, tf_offsets.data(), sz_tf_off);
        memcpy(hp + off_tf_min, tf_mins.data(), sz_tf_min);
        memcpy(hp + off_tf_cf, tf_cont_fs.data(), sz_tf_cf);
        CUDA_CHECK(cudaMemcpyAsync(static_cast<char*>(impl->all_d_tf_pack.ptr), hp, tf_pack_sz, cudaMemcpyHostToDevice, stream));
    }

    float* h_sc_out = nullptr;
    if (n_sc > 0) {
        // sc_disc is already pinned (sideChainAnglesDiscretized_flat uses CudaPinnedAllocator)
        if (!out_d_sc) h_sc_out = impl->all_ph_sc_out.ensure_n<float>(n_sc);
        CUDA_CHECK(cudaMemcpyAsync(d_sc, sc_disc.data(), n_sc * sizeof(uint8_t), cudaMemcpyHostToDevice, stream));
    }

    // --- Backbone per-struct params: pack [offsets | mins | cont_fs] into one device buffer ---
    // bb_struct_offsets: n_structs+1 ints; bb_mins_all, bb_cont_fs_all: n_structs*6 floats each
    const size_t sz_bb_off = (size_t)(n_structs + 1) * sizeof(int);
    const size_t sz_bb_min = (size_t)n_structs * 6 * sizeof(float);
    const size_t sz_bb_cf = (size_t)n_structs * 6 * sizeof(float);
    const size_t off_bb_min = sz_bb_off;
    const size_t off_bb_cf = off_bb_min + sz_bb_min;
    const size_t bb_params_sz = off_bb_cf + sz_bb_cf;

    char* dp_bb = static_cast<char*>(impl->all_d_bb_params.ensure(bb_params_sz, stream));
    char* hp_bb = static_cast<char*>(impl->all_ph_bb_params.ensure(bb_params_sz));
    memcpy(hp_bb, bb_struct_offsets.data(), sz_bb_off);
    memcpy(hp_bb + off_bb_min, bb_mins_all.data(), sz_bb_min);
    memcpy(hp_bb + off_bb_cf, bb_cont_fs_all.data(), sz_bb_cf);
    CUDA_CHECK(cudaMemcpyAsync(dp_bb, hp_bb, bb_params_sz, cudaMemcpyHostToDevice, stream));

    auto* d_bb_offsets = reinterpret_cast<int*>(dp_bb);
    auto* d_bb_mins = reinterpret_cast<float*>(dp_bb + off_bb_min);
    auto* d_bb_cf = reinterpret_cast<float*>(dp_bb + off_bb_cf);

    int blk_bb = (n_bb_residues + 255) / 256;
    unpack_raw_file_continuize_bb_per_struct_kernel<<<blk_bb, 256, 0, stream>>>(
        d_bb_raw, d_bb_out, d_bb_offsets, d_bb_mins, d_bb_cf, n_structs, n_bb_residues);

    // --- Launch TF kernel (uint8_t input, per-struct calibration) ---
    if (n_tf > 0 && n_structs > 0) {
        int blk_tf = (n_tf + 255) / 256;
        continuize_per_struct_u8_kernel<<<blk_tf, 256, 0, stream>>>(
            d_tf, d_tf_out, d_tf_off, d_tf_min, d_tf_cf, n_structs, n_tf);
    }

    // --- Launch SC kernel (uint8_t input, single global calibration) ---
    if (n_sc > 0) {
        int blk_sc = (n_sc + 255) / 256;
        continuize_fused_u8_kernel<<<blk_sc, 256, 0, stream>>>(d_sc, d_sc_out, n_sc, sc_min, sc_cont_f);
    }

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
        throw std::runtime_error("continuize_all kernel failed: " + std::string(cudaGetErrorString(err)));

    // --- Download results ---
    // bb_cont is pinned (PinnedFloatVec), so cudaMemcpyAsync can DMA directly to it
    if (!out_d_bb) {
        bb_cont.resize(n_bb_out);
        CUDA_CHECK(cudaMemcpyAsync(bb_cont.data(), d_bb_out, n_bb_out * sizeof(float), cudaMemcpyDeviceToHost, stream));
    }
    if (h_tf_out) CUDA_CHECK(cudaMemcpyAsync(h_tf_out, d_tf_out, n_tf * sizeof(float), cudaMemcpyDeviceToHost, stream));
    if (h_sc_out) CUDA_CHECK(cudaMemcpyAsync(h_sc_out, d_sc_out, n_sc * sizeof(float), cudaMemcpyDeviceToHost, stream));

    // No sync here: all D→H transfers above are enqueued on the stream and will
    // complete before the caller's final cudaEventSynchronize (inside fillAtomCoordinatesGPU).
    // FIFO stream ordering guarantees h_tf_out is valid by the time OXT fill reads it.

    if (out_d_bb) {
        *out_d_bb = d_bb_out;
        bb_cont.clear();
    }
    if (out_h_tf) *out_h_tf = h_tf_out;
    if (out_d_sc) *out_d_sc = d_sc_out;
    if (out_d_tf) *out_d_tf = d_tf_out;
}
