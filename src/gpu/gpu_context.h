// SPDX-License-Identifier: MIT
/**
 * File: gpu_context.h
 * Project: foldcomp
 * Description:
 *     Persistent GPU buffer context structs for backbone and sidechain modules.
 *     Replaces static/global storage with explicit per-Foldcomp-object lifetime.
 */
#pragma once
#include "gpu_buffer.h"
#include "gpu_stream.h"
#include <cstddef>
#include <cstdint>

/**
 * All persistent device/pinned buffers used by the backbone reconstruction
 * functions in gpu_nerf.cu. One instance per Foldcomp object.
 */
struct NerfBuffers {
    // Cached device pointers set by reconstructBackboneGPU_segmented,
    // consumed by buildBackboneForSidechains on the next call.
    float* last_bb_dev_ptr = nullptr;
    const int8_t* last_d_rt = nullptr;

    // buildBackboneForSidechains: single pack buffer [pa | rfo | bbo]
    GrowOnlyDeviceBuf d_bb_sc_pack;
    GrowOnlyPinnedBuf ph_bb_sc_pack;

    // reconstructBackboneGPU_segmented working buffers
    GrowOnlyDeviceBuf unit_prev_atoms; // [total_units * 9] per-unit sliding-window state
    GrowOnlyDeviceBuf output;          // [total_output floats] backbone coordinate output
    GrowOnlyDeviceBuf seg_pack;        // packed per-protein metadata (single H→D transfer)
    GrowOnlyPinnedBuf ph_seg_pack;     // pinned staging for seg_pack

    /**
     * Pre-allocate pinned buffers to their expected steady-state sizes.
     * Call once after CUDA init, before the first batch, to avoid cold-start
     * cudaMallocHost latency on the first few batches.
     *
     * @param max_batch          Maximum number of structures per batch
     * @param max_total_residues Upper bound on total residues across the batch
     *                           (max_batch * max_residues_per_struct)
     */
    /**
     * @param max_batch          Maximum number of structures per batch
     * @param max_total_residues Upper bound on total residues across the batch
     */
    void preallocate(int max_batch, int max_total_residues, cudaStream_t stream = 0) {
        const size_t sz_bb_sc_pack =
            (size_t)max_batch * 9 * sizeof(float) +
            (size_t)(max_batch + 1) * 2 * sizeof(int);

        // anchor count ≈ max_total_residues / 25 (DEFAULT_ANCHOR_THRESHOLD)
        const size_t max_anchors = (size_t)max_total_residues / 25 + max_batch;
        const size_t sz_seg_pack =
            (size_t)max_batch * 9 * sizeof(float) +          // prev_atoms
            (size_t)(max_batch + 1) * sizeof(int) +          // aoff
            (size_t)(max_total_residues + 1) * sizeof(int) + // roff
            max_anchors * sizeof(int) +                      // aidx
            (size_t)(max_batch + 1) * sizeof(int) +          // aoffg
            max_anchors * 9 * sizeof(float) +                // acoords
            (size_t)(max_batch + 1) * sizeof(int) +          // acoff
            (size_t)max_batch * sizeof(int) +                // nanch
            (size_t)max_batch * sizeof(int) +                // nres
            (size_t)(max_batch + 1) * sizeof(int) +          // coff
            (size_t)max_total_residues * sizeof(int8_t);     // rcodes

        // Pinned host buffers
        ph_bb_sc_pack.preallocate(sz_bb_sc_pack);
        ph_seg_pack.preallocate(sz_seg_pack);

        // Mirror device buffers (avoid batch-N cudaMallocAsync stalls mid-run)
        d_bb_sc_pack.preallocate(sz_bb_sc_pack, stream);
        seg_pack.preallocate(sz_seg_pack, stream);
        // output: sum of (n_residues-1)*9 floats per protein ≈ max_total_residues * 9
        output.preallocate((size_t)max_total_residues * 9 * sizeof(float), stream);
        // unit_prev_atoms: n_proteins * max_segments * 9; max_segments ≈ max_anchors / max_batch
        unit_prev_atoms.preallocate(max_anchors * 9 * sizeof(float), stream);
    }
};

/**
 * All persistent device/pinned buffers used by sidechain reconstruction
 * functions in gpu_sidechain.cu. One instance per Foldcomp object.
 */
struct SidechainBuffers {
    // reconstructSidechainsGPU output buffers (consumed by fillAtomCoordinatesGPU)
    GrowOnlyDeviceBuf output_coords;     // [total_atoms_ub * 3 floats]
    GrowOnlyDeviceBuf output_atom_names; // [total_atoms_ub int8_t]

    // Cached device pointers set by reconstructSidechainsGPU,
    // read by fillAtomCoordinatesGPU (same stream, no sync needed).
    const int8_t* d_rt_cached = nullptr;
    const int2* d_ooff_cached = nullptr; // .x=sc angle offset, .y=atom output offset

    // buildSidechainMeta: prefix-sum input/output buffers (int2: .x=sc count/off, .y=atom count/off)
    GrowOnlyDeviceBuf d_meta_counts; // int2[n_residues+1]
    GrowOnlyDeviceBuf d_meta_off;    // int2[n_residues+1]
    GrowOnlyDeviceBuf d_cub_temp;    // CUB temporary storage

    // buildSidechainMeta: bucket classification buffers
    GrowOnlyDeviceBuf d_bucket_indices; // int[n_residues] — all 4 buckets contiguous
    GrowOnlyDeviceBuf d_bsz_boff_bcnt;  // int[12]: bsz[4] | boff[4] | bcnt[4]

    // reconstructSidechainsGPU: use_alt_per_struct upload (≤100 bytes)
    GrowOnlyDeviceBuf d_alt;
    GrowOnlyPinnedBuf ph_alt;

    // fillAtomCoordinatesGPU: per-struct metadata pack + AtomCoordinate output
    GrowOnlyDeviceBuf d_fill_pack;
    GrowOnlyPinnedBuf ph_fill_pack;
    GrowOnlyDeviceBuf d_ac_output; // AtomCoordinate[total_output_atoms]

    /**
     * @param max_batch          Maximum number of structures per batch
     * @param max_total_residues Upper bound on total residues across the batch
     * @param max_total_atoms    Upper bound on total heavy atoms across the batch
     *                           (max_total_atoms_bytes passed separately for d_ac_output
     *                            since sizeof(AtomCoordinate) is not known here)
     * @param ac_output_bytes    sizeof(AtomCoordinate) * max_total_atoms
     */
    void preallocate(int max_batch, int max_total_residues, int max_total_atoms,
        size_t ac_output_bytes, cudaStream_t stream = 0) {
        const size_t sz_fill_pack =
            (size_t)(max_batch + 1) * sizeof(int) +
            (size_t)max_batch * sizeof(int) +
            (size_t)max_batch * sizeof(int) +
            (size_t)max_batch * sizeof(char);

        // Pinned host buffers
        ph_alt.preallocate((size_t)max_batch * sizeof(int8_t));
        ph_fill_pack.preallocate(sz_fill_pack);

        // Device buffers
        d_alt.preallocate((size_t)max_batch * sizeof(int8_t), stream);
        d_fill_pack.preallocate(sz_fill_pack, stream);
        output_coords.preallocate((size_t)max_total_atoms * 3 * sizeof(float), stream);
        output_atom_names.preallocate((size_t)max_total_atoms * sizeof(int8_t), stream);
        d_meta_counts.preallocate((size_t)(max_total_residues + 1) * sizeof(int2), stream);
        d_meta_off.preallocate((size_t)(max_total_residues + 1) * sizeof(int2), stream);
        d_bucket_indices.preallocate((size_t)max_total_residues * sizeof(int), stream);
        d_ac_output.preallocate(ac_output_bytes, stream);
    }
};
