// SPDX-License-Identifier: MIT

/**
 * File: gpu_sidechain.h
 * Project: foldcomp
 * Created: 2026-01-31
 * Description:
 *     GPU-accelerated sidechain reconstruction using CUDA.
 *     Provides batch processing of amino acid sidechain atoms using NeRF algorithm.
 */

#pragma once

#include "atom_coordinate.h"
#include "gpu_context.h"
#include <cstdint>
#include <vector>

// ============================================================================
// GPU-friendly Amino Acid Data Structures
// ============================================================================

// Maximum atoms per amino acid (TRP has 14 atoms)
static constexpr int GPU_MAX_ATOMS_PER_AA = 14;

// Maximum sidechain atoms (total - 3 backbone atoms)
static constexpr int GPU_MAX_SIDECHAIN_ATOMS = 11;

// Number of standard amino acids
static constexpr int GPU_NUM_AMINO_ACIDS = 20;

// Atom name indices (for encoding atom names as integers)
enum AtomNameIndex : int8_t {
    ATOM_N = 0,
    ATOM_CA = 1,
    ATOM_C = 2,
    ATOM_O = 3,
    ATOM_CB = 4,
    ATOM_CG = 5,
    ATOM_CG1 = 6,
    ATOM_CG2 = 7,
    ATOM_CD = 8,
    ATOM_CD1 = 9,
    ATOM_CD2 = 10,
    ATOM_CE = 11,
    ATOM_CE1 = 12,
    ATOM_CE2 = 13,
    ATOM_CE3 = 14,
    ATOM_CZ = 15,
    ATOM_CZ2 = 16,
    ATOM_CZ3 = 17,
    ATOM_CH2 = 18,
    ATOM_ND1 = 19,
    ATOM_ND2 = 20,
    ATOM_NE = 21,
    ATOM_NE1 = 22,
    ATOM_NE2 = 23,
    ATOM_NH1 = 24,
    ATOM_NH2 = 25,
    ATOM_NZ = 26,
    ATOM_OD1 = 27,
    ATOM_OD2 = 28,
    ATOM_OE1 = 29,
    ATOM_OE2 = 30,
    ATOM_OG = 31,
    ATOM_OG1 = 32,
    ATOM_OH = 33,
    ATOM_SD = 34,
    ATOM_SG = 35,
    ATOM_INVALID = -1
};

// Amino acid type indices (matching AminoAcidIndexMap order)
enum AminoAcidIndex : int8_t {
    AA_ALA = 0,
    AA_ARG = 1,
    AA_ASN = 2,
    AA_ASP = 3,
    AA_CYS = 4,
    AA_GLN = 5,
    AA_GLU = 6,
    AA_GLY = 7,
    AA_HIS = 8,
    AA_ILE = 9,
    AA_LEU = 10,
    AA_LYS = 11,
    AA_MET = 12,
    AA_PHE = 13,
    AA_PRO = 14,
    AA_SER = 15,
    AA_THR = 16,
    AA_TRP = 17,
    AA_TYR = 18,
    AA_VAL = 19,
    AA_UNK = 20
};

/**
 * @brief GPU-friendly sidechain atom description
 *
 * Describes how to build a sidechain atom using NeRF algorithm.
 * Each atom depends on 3 previous atoms.
 */
struct GPUSidechainAtom {
    int8_t atom_name;    // AtomNameIndex for output atom
    int8_t dep_atoms[3]; // AtomNameIndex for 3 dependency atoms
    float bond_length;   // Bond length from dep_atoms[2] to this atom
    float bond_angle;    // Angle at dep_atoms[1]-dep_atoms[2]-this (degrees)
};

/**
 * @brief GPU-friendly amino acid structure
 *
 * Contains all information needed to reconstruct an amino acid on GPU.
 * Uses integer indices instead of strings for GPU compatibility.
 */
struct GPUAminoAcid {
    int8_t n_atoms;                                      // Total number of atoms
    int8_t n_sidechain_atoms;                            // Number of sidechain atoms (n_atoms - 3)
    int8_t atom_order[GPU_MAX_ATOMS_PER_AA];             // Atom names in order (AtomNameIndex)
    int8_t alt_atom_order[GPU_MAX_ATOMS_PER_AA];         // Alternate atom order for reordering
    GPUSidechainAtom sidechain[GPU_MAX_SIDECHAIN_ATOMS]; // Sidechain atom build instructions
};

// ============================================================================
// Host-side Functions
// ============================================================================

/**
 * @brief Initialize GPU amino acid lookup table
 *
 * Must be called once per device before using any GPU sidechain functions
 * on that device. Copies amino acid data to the current device's constant
 * memory; tracks initialization per-device internally since constant memory
 * is not shared across devices.
 */
void initGPUAminoAcidTable();

/**
 * @brief Clean up GPU amino acid resources
 *
 * Resets initialization tracking for all devices.
 */
void cleanupGPUAminoAcidTable();

/**
 * @brief Check if GPU amino acid table is initialized on the current device
 */
bool isGPUAminoAcidTableInitialized();

/**
 * @brief Convert 3-letter amino acid code to index
 * @param code Three-letter code (e.g., "ALA", "ARG")
 * @return AminoAcidIndex value, or AA_UNK if unknown
 */
AminoAcidIndex getAminoAcidIndex(const char* code);

/**
 * @brief Convert single-letter code to amino acid index
 * @param code Single letter code (e.g., 'A', 'R')
 * @return AminoAcidIndex value, or AA_UNK if unknown
 */
AminoAcidIndex getAminoAcidIndexFromChar(char code);

/**
 * @brief Get sidechain torsion count directly from BackboneChain.residue int (5-bit encoding)
 * @param residue_int BackboneChain.residue value (0–19 for standard AAs)
 * @return Torsion count, or 0 for unknown
 */
int getSCTorsionCountFromInt(unsigned int residue_int);

/**
 * @brief Get total atom count directly from BackboneChain.residue int (5-bit encoding)
 * @param residue_int BackboneChain.residue value (0–19 for standard AAs)
 * @return Total atom count (including backbone), or 4 (GLY) for unknown
 */
int getAtomCountFromInt(unsigned int residue_int);

/**
 * @brief Get number of sidechain torsion angles for an amino acid
 * @param aa_idx Amino acid index
 * @return Number of torsion angles needed
 */
int getGPUSidechainTorsionCount(AminoAcidIndex aa_idx);

/**
 * @brief Get number of atoms for an amino acid
 * @param aa_idx Amino acid index
 * @return Total number of atoms
 */
int getGPUAtomCount(AminoAcidIndex aa_idx);

/**
 * @brief Convert AtomNameIndex to string
 * @param idx Atom name index
 * @return Atom name string (e.g., "CA", "CB")
 */
const char* getAtomNameString(AtomNameIndex idx);

// ============================================================================
// GPU Sidechain Meta (computed from on-device residue types)
// ============================================================================

/**
 * @brief Device pointers to sidechain metadata arrays, computed entirely on GPU
 *        from the on-device residue types via buildSidechainMeta().
 *
 * d_sc_angle_offsets and d_output_offsets are exclusive prefix sums of per-residue
 * SC torsion counts and atom counts, respectively.  Both have n_residues+1 elements.
 * d_bsz/d_boff/d_bidx hold per-bucket metadata entirely on device. The templated
 * sidechain kernels read bucket sizes and offsets directly from device memory.
 */
struct SCMetaPtrs {
    const int2* d_offsets; // [n_residues+1] .x=sc angle offset, .y=atom output offset
    const int* d_bsz;      // device int[4]: residue count per bucket
    const int* d_boff;     // device int[4]: start offset of each bucket in d_bidx
    const int* d_bidx;     // device int[n_residues]: scatter index array (all 4 buckets)
};

/**
 * @brief Build sidechain metadata on GPU from already-on-device residue types.
 *
 * Computes sc_angle_offsets and output_offsets via CUB exclusive prefix sums,
 * and always classifies residues into 4 buckets using GPU atomics.
 *
 * @param d_rt         Device pointer to residue types [n_residues] (from BBSidechainPtrs.d_rt)
 * @param n_residues   Total residues across the batch
 * @return SCMetaPtrs with all device pointers valid until the next buildSidechainMeta call
 */
SCMetaPtrs buildSidechainMeta(SidechainBuffers& ctx, const int8_t* d_rt, int n_residues, cudaStream_t stream);

// ============================================================================
// GPU Batch Sidechain Reconstruction
// ============================================================================

/**
 * @brief Reconstruct sidechains for a batch of residues on GPU
 *
 * This is the main entry point for GPU sidechain reconstruction.
 * Processes all residues in parallel on the GPU.
 *
 * @param n_residues Total number of residues to process
 * @param residue_types Array of amino acid indices [n_residues]
 * @param backbone_coords Backbone atom coordinates (N, CA, C per residue)
 *                        Layout: [res0_N_x, res0_N_y, res0_N_z, res0_CA_x, ..., res0_C_z, res1_N_x, ...]
 *                        Size: n_residues * 9 floats
 * @param sidechain_angles Flattened sidechain torsion angles
 *                         Variable length per residue based on residue type
 * @param sc_angle_offsets Offset into sidechain_angles for each residue [n_residues + 1]
 * @param use_alt_order Whether to use alternate atom ordering [n_residues]
 * @param output_coords Output array for all atom coordinates
 *                      Layout: [res0_atom0_x, y, z, res0_atom1_x, y, z, ...]
 * @param output_offsets Offset into output_coords for each residue [n_residues + 1]
 * @param output_atom_names Output atom name indices (for building AtomCoordinate objects)
 *                          Size: total atoms across all residues
 */
/**
 * Backbone is supplied as four already-on-device pointers from buildBackboneForSidechains.
 * Residue types, sc_angle_offsets, and output_offsets are supplied via device pointers in
 * SCMetaPtrs (from buildSidechainMeta — no CPU→GPU upload needed for those arrays).
 *
 * use_alt_per_struct (one byte per struct, max 100 bytes) is the only remaining H→D upload.
 * d_sidechain_angles must be the device pointer returned by continuize_all.
 */
void reconstructSidechainsGPU(
    SidechainBuffers& ctx,
    int n_residues,
    int n_structs,
    const int8_t* d_rt,
    const float* d_bb_out,
    const float* d_pa,
    const int* d_rfo,
    const int* d_bbo,
    const float* d_sidechain_angles,
    const std::vector<int8_t>& use_alt_per_struct,
    const SCMetaPtrs& meta,
    cudaStream_t stream);

/**
 * @brief Fill an AtomCoordinate array on GPU using on-device sidechain outputs.
 *
 * Must be called immediately after reconstructSidechainsGPU (outputs remain on device).
 * Uploads small per-struct metadata arrays, runs a per-residue fill kernel, downloads result.
 * Per-residue values (local_idx, atom_num_start, output_base, chain) are computed in the kernel
 * from per-struct data, eliminating O(total_residues) CPU prep and ~700 KB of H→D uploads.
 *
 * @param n_residues            Total residues across the batch
 * @param n_structs             Number of structures in the batch
 * @param total_atoms           Total atoms (= output_offsets.back() from reconstructSidechainsGPU)
 * @param total_output_atoms    >= total_atoms (includes OXT slots)
 * @param n_tf                  Length of d_tf (typically == n_residues)
 * @param d_tf                  Device ptr to continuized TF floats [n_tf] (from continuize_all)
 * @param struct_res_offsets    Cumulative residue counts [n_structs + 1]
 * @param base_anum_per_struct  header.idxAtom per structure [n_structs]
 * @param base_ridx_per_struct  header.idxResidue per structure [n_structs]
 * @param base_delta_per_struct ac_offset - gpu_off_start per structure [n_structs]
 * @param chain_per_struct      Chain bytes per structure [n_structs * CHAIN_ID_LENGTH], mirrors ChainId (FixedStr<CHAIN_ID_LENGTH>)
 * @param output                Pre-allocated AtomCoordinate array [total_output_atoms]
 */
/**
 * This function no longer blocks.  After recording e_fill_done on stream A,
 * it launches the D→H on stream B and records e_slot_ready and cpu_sync_event.
 * The caller must cudaEventSynchronize(cpu_sync_event) before reading `output`.
 */
void fillAtomCoordinatesGPU(
    SidechainBuffers& ctx,
    int n_residues,
    int n_structs,
    int total_atoms,
    int total_output_atoms,
    int n_tf,
    const float* d_tf,
    const std::vector<int>& struct_res_offsets,
    const std::vector<int>& base_anum_per_struct,
    const std::vector<int>& base_ridx_per_struct,
    const std::vector<int>& base_delta_per_struct,
    const std::vector<char>& chain_per_struct,
    AtomCoordinate* output,
    cudaStream_t stream_a,     // compute stream for this pipeline
    cudaStream_t stream_b,     // D→H stream for this pipeline
    cudaEvent_t e_fill_done,   // recorded after fill kernel, signals stream_b
    cudaEvent_t e_slot_ready,  // recorded after D→H, guards next reuse of output buffer
    cudaEvent_t cpu_sync_event,// per-batch: recorded after D→H for CPU sync
    bool skip_download = false // when true, run the fill kernel but skip the D→H
                               // (caller consumes ctx.d_ac_output on the device)
);

//#endif // GPU_SIDECHAIN_H is already defined by pragma once
