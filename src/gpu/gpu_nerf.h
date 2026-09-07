// SPDX-License-Identifier: MIT

/**
 * File: gpu_nerf.h
 * Project: foldcomp
 * Created: 2026-01-06
 * Description:
 *     GPU-accelerated NeRF (Natural Extension of Reference Frame) operations using CUDA.
 *     Implements anchor-segment parallel backbone reconstruction.
 */

#pragma once
#include "gpu_buffer.h"
#include "gpu_context.h"
#include <vector>
#include <cstdint>
#include <cuda_runtime.h>

/**
 * @brief Reconstruct backbone atom coordinates for a batch of structures via
 * anchor-segment parallel NeRF placement.
 *
 * Each structure's residues are split into independent segments anchored at
 * known coordinates (anchor_indices_flat/anchor_coords_flat), allowing all
 * segments across the batch to be placed in parallel on the GPU.
 *
 * @param ctx                       GPU buffer context (device/pinned scratch storage)
 * @param n_proteins                Number of structures in the batch
 * @param n_residues_per_protein    Residue count per structure [n_proteins]
 * @param prev_atoms_flat           First 3 atoms per structure [n_proteins * 9]
 * @param backbone_offsets          Cumulative backbone-chain offsets [n_proteins+1]
 * @param residue_types             Residue type per residue, flattened across the batch
 * @param residue_offsets           Cumulative residue offsets [n_proteins+1]
 * @param anchor_indices_flat       Anchor residue indices, flattened across the batch
 * @param anchor_offsets            Cumulative anchor-index offsets [n_proteins+1]
 * @param anchor_coords_flat        Anchor xyz coordinates, flattened across the batch
 * @param anchor_coord_offsets      Cumulative anchor-coordinate offsets [n_proteins+1]
 * @param n_anchors_per_protein     Anchor count per structure [n_proteins]
 * @param output_coords             Output: reconstructed xyz coordinates (pinned host)
 * @param coord_offsets             Output: cumulative output-coordinate offsets [n_proteins+1]
 * @param d_angles                  Device pointer to continuized backbone torsion angles
 * @param stream                    CUDA stream to launch on
 */
void reconstructBackboneGPU_segmented(
    NerfBuffers& ctx,
    int n_proteins,
    const std::vector<int>& n_residues_per_protein,
    const std::vector<float>& prev_atoms_flat,
    const std::vector<int>& backbone_offsets,
    const std::vector<int8_t>& residue_types,
    const std::vector<int>& residue_offsets,
    const std::vector<int>& anchor_indices_flat,
    const std::vector<int>& anchor_offsets,
    const std::vector<float>& anchor_coords_flat,
    const std::vector<int>& anchor_coord_offsets,
    const std::vector<int>& n_anchors_per_protein,
    PinnedFloatVec& output_coords,
    std::vector<int>& coord_offsets,
    const float* d_angles,
    cudaStream_t stream);

/**
 * @brief Device pointers returned by buildBackboneForSidechains.
 *
 * All pointers remain valid until the next buildBackboneForSidechains call.
 * Pass directly to reconstructSidechainsGPU — no additional H→D transfer needed.
 */
struct BBSidechainPtrs {
    float* d_bb_out;    // backbone recon output (residues 1..n-1 per struct)
    float* d_pa;        // prevAtoms per struct [n_structs * 9]
    int* d_rfo;         // residue full offsets [n_structs+1]: prefix sums of n_residues
    int* d_bbo;         // bb output float offsets [n_structs+1]
    const int8_t* d_rt; // residue types [total_residues] — pointer into backbone static buf
};

/**
 * @brief Upload backbone metadata for sidechain reconstruction.
 *
 * Must be called immediately after reconstructBackboneGPU_segmented.
 * Uploads a single pack [prevAtoms | rfo | bbo] to device and returns device pointers.
 * No scatter kernel is launched; backbone lookup is inlined into the sidechain kernels.
 *
 * @param n_structs             Number of structures in batch
 * @param coord_offsets_cpu     Backbone float offsets [n_structs+1] (from backbone recon)
 * @param n_residues_per_protein Residue counts [n_structs]
 * @param prev_atoms_flat       First 3 atoms per structure [n_structs * 9]
 * @param total_residues        (unused, kept for caller compat)
 * @return BBSidechainPtrs with all device pointers valid until next call
 */
BBSidechainPtrs buildBackboneForSidechains(
    NerfBuffers& ctx,
    int n_structs,
    const std::vector<int>& coord_offsets_cpu,
    const std::vector<int>& n_residues_per_protein,
    const std::vector<float>& prev_atoms_flat,
    int total_residues,
    cudaStream_t stream);
