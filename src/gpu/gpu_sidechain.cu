// SPDX-License-Identifier: MIT

/**
 * File: gpu_sidechain.cu
 * Project: foldcomp
 * Created: 2026-01-31
 * Description:
 *     GPU-accelerated sidechain reconstruction using CUDA.
 *     Implements NeRF algorithm for parallel sidechain atom placement.
 */

#include "gpu_sidechain.h"
#include "gpu_context.h"
#include "gpu_buffer.h"
#include "gpu_stream.h"
#include "atom_coordinate.h"
#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <cmath>
#include <nvtx3/nvtx3.hpp>
#include <cstring>
#include <stdexcept>
#include <iostream>
#include <mutex>
#include <unordered_set>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

// ============================================================================
// GPU Amino Acid Lookup Table (Constant Memory)
// ============================================================================

// Store amino acid data in constant memory for fast access
__constant__ GPUAminoAcid d_amino_acids[GPU_NUM_AMINO_ACIDS + 1]; // +1 for UNK

// Constant memory (d_amino_acids, d_atom_name_strs, etc.) is per-device, so
// initialization must be tracked per-device rather than with a single flag.
static std::mutex g_aa_table_mutex;
static std::unordered_set<int> g_aa_table_initialized_devices;

// ============================================================================
// Constant memory string lookup tables for AtomCoordinate fill kernel
// ============================================================================
// Atom name strings: 36 entries (ATOM_N..ATOM_SG), max 4 chars + null
__constant__ char d_atom_name_strs[36][5];
// Residue name strings: 21 entries (AA_ALA..AA_UNK), max 3 chars + null
__constant__ char d_residue_name_strs[21][4];

// Per-AA lookup tables for buildSidechainMeta kernels (indexed by AminoAcidIndex 0..20)
// int8_t: max value = 11 (TRP torsion count) — fits in signed byte
__constant__ int8_t d_sc_torsion_counts_k[21]; // SIDECHAIN_TORSION_COUNTS per AA
__constant__ int8_t d_atom_counts_k[21];       // ATOM_COUNTS per AA
__constant__ int8_t d_bucket_id_k[21];         // AA_BUCKET_ID per AA (0-3)

// ============================================================================
// Host-side Amino Acid Data Initialization
// ============================================================================

// Atom name lookup table
static const char* ATOM_NAME_STRINGS[] = {
    "N", "CA", "C", "O", "CB", "CG", "CG1", "CG2", "CD", "CD1", "CD2",
    "CE", "CE1", "CE2", "CE3", "CZ", "CZ2", "CZ3", "CH2",
    "ND1", "ND2", "NE", "NE1", "NE2", "NH1", "NH2", "NZ",
    "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "SD", "SG"};

const char* getAtomNameString(AtomNameIndex idx) {
    if (idx < 0 || idx > ATOM_SG) return "UNK";
    return ATOM_NAME_STRINGS[idx];
}

// Sidechain torsion counts per amino acid (from getSideChainTorsionNum)
static const int SIDECHAIN_TORSION_COUNTS[GPU_NUM_AMINO_ACIDS] = {
    2,  // ALA
    8,  // ARG
    5,  // ASN
    5,  // ASP
    3,  // CYS
    6,  // GLN
    6,  // GLU
    1,  // GLY
    7,  // HIS
    5,  // ILE
    5,  // LEU
    6,  // LYS
    5,  // MET
    8,  // PHE
    4,  // PRO
    3,  // SER
    4,  // THR
    11, // TRP
    9,  // TYR
    4   // VAL
};

int getSCTorsionCountFromInt(unsigned int residue_int) {
    if (residue_int >= (unsigned int)GPU_NUM_AMINO_ACIDS) return 0;
    return SIDECHAIN_TORSION_COUNTS[residue_int];
}

// Atom counts per amino acid
static const int ATOM_COUNTS[GPU_NUM_AMINO_ACIDS] = {
    5,  // ALA: N, CA, C, O, CB
    11, // ARG: N, CA, C, O, CB, CG, CD, NE, CZ, NH1, NH2
    8,  // ASN: N, CA, C, O, CB, CG, OD1, ND2
    8,  // ASP: N, CA, C, O, CB, CG, OD1, OD2
    6,  // CYS: N, CA, C, O, CB, SG
    9,  // GLN: N, CA, C, O, CB, CG, CD, OE1, NE2
    9,  // GLU: N, CA, C, O, CB, CG, CD, OE1, OE2
    4,  // GLY: N, CA, C, O
    10, // HIS: N, CA, C, O, CB, CG, ND1, CD2, CE1, NE2
    8,  // ILE: N, CA, C, O, CB, CG1, CG2, CD1
    8,  // LEU: N, CA, C, O, CB, CG, CD1, CD2
    9,  // LYS: N, CA, C, O, CB, CG, CD, CE, NZ
    8,  // MET: N, CA, C, O, CB, CG, SD, CE
    11, // PHE: N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ
    7,  // PRO: N, CA, C, O, CB, CG, CD
    6,  // SER: N, CA, C, O, CB, OG
    7,  // THR: N, CA, C, O, CB, OG1, CG2
    14, // TRP: N, CA, C, O, CB, CG, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3, CH2
    12, // TYR: N, CA, C, O, CB, CG, CD1, CD2, CE1, CE2, CZ, OH
    7   // VAL: N, CA, C, O, CB, CG1, CG2
};

int getGPUSidechainTorsionCount(AminoAcidIndex aa_idx) {
    if (aa_idx < 0 || aa_idx >= GPU_NUM_AMINO_ACIDS) return 0;
    return SIDECHAIN_TORSION_COUNTS[aa_idx];
}

int getGPUAtomCount(AminoAcidIndex aa_idx) {
    // AA_UNK, and anything else out of range, is UNK: N, CA, C only — see
    // h_amino_acids[AA_UNK] in buildHostAminoAcidTable.
    if (aa_idx < 0 || aa_idx >= GPU_NUM_AMINO_ACIDS) return 3;
    return ATOM_COUNTS[aa_idx];
}

int getAtomCountFromInt(unsigned int residue_int) {
    // Any raw on-disk residue code outside the known table (including the
    // encoder's own "unknown residue" sentinel, AA_UNK_INT = 23 — see
    // utility.h — which does not equal GPU's AA_UNK = 20) is UNK: N, CA, C
    // only. See h_amino_acids[AA_UNK] in buildHostAminoAcidTable.
    if (residue_int >= (unsigned int)GPU_NUM_AMINO_ACIDS) return 3;
    return ATOM_COUNTS[residue_int];
}

AminoAcidIndex getAminoAcidIndex(const char* code) {
    if (!code) return AA_UNK;
    if (strcmp(code, "ALA") == 0) return AA_ALA;
    if (strcmp(code, "ARG") == 0) return AA_ARG;
    if (strcmp(code, "ASN") == 0) return AA_ASN;
    if (strcmp(code, "ASP") == 0) return AA_ASP;
    if (strcmp(code, "CYS") == 0) return AA_CYS;
    if (strcmp(code, "GLN") == 0) return AA_GLN;
    if (strcmp(code, "GLU") == 0) return AA_GLU;
    if (strcmp(code, "GLY") == 0) return AA_GLY;
    if (strcmp(code, "HIS") == 0) return AA_HIS;
    if (strcmp(code, "ILE") == 0) return AA_ILE;
    if (strcmp(code, "LEU") == 0) return AA_LEU;
    if (strcmp(code, "LYS") == 0) return AA_LYS;
    if (strcmp(code, "MET") == 0) return AA_MET;
    if (strcmp(code, "PHE") == 0) return AA_PHE;
    if (strcmp(code, "PRO") == 0) return AA_PRO;
    if (strcmp(code, "SER") == 0) return AA_SER;
    if (strcmp(code, "THR") == 0) return AA_THR;
    if (strcmp(code, "TRP") == 0) return AA_TRP;
    if (strcmp(code, "TYR") == 0) return AA_TYR;
    if (strcmp(code, "VAL") == 0) return AA_VAL;
    return AA_UNK;
}

AminoAcidIndex getAminoAcidIndexFromChar(char code) {
    switch (code) {
    case 'A':
        return AA_ALA;
    case 'R':
        return AA_ARG;
    case 'N':
        return AA_ASN;
    case 'D':
        return AA_ASP;
    case 'C':
        return AA_CYS;
    case 'Q':
        return AA_GLN;
    case 'E':
        return AA_GLU;
    case 'G':
        return AA_GLY;
    case 'H':
        return AA_HIS;
    case 'I':
        return AA_ILE;
    case 'L':
        return AA_LEU;
    case 'K':
        return AA_LYS;
    case 'M':
        return AA_MET;
    case 'F':
        return AA_PHE;
    case 'P':
        return AA_PRO;
    case 'S':
        return AA_SER;
    case 'T':
        return AA_THR;
    case 'W':
        return AA_TRP;
    case 'Y':
        return AA_TYR;
    case 'V':
        return AA_VAL;
    default:
        return AA_UNK;
    }
}

// Build amino acid table on host
static void buildHostAminoAcidTable(GPUAminoAcid* h_amino_acids) {
    memset(h_amino_acids, 0, sizeof(GPUAminoAcid) * (GPU_NUM_AMINO_ACIDS + 1));

// Helper macro for cleaner initialization
#define SET_SC_ATOM(aa, idx, name, d0, d1, d2, bl, ba) \
    h_amino_acids[aa].sidechain[idx].atom_name = name; \
    h_amino_acids[aa].sidechain[idx].dep_atoms[0] = d0; \
    h_amino_acids[aa].sidechain[idx].dep_atoms[1] = d1; \
    h_amino_acids[aa].sidechain[idx].dep_atoms[2] = d2; \
    h_amino_acids[aa].sidechain[idx].bond_length = bl; \
    h_amino_acids[aa].sidechain[idx].bond_angle = ba;

    // ALA - Alanine
    h_amino_acids[AA_ALA].n_atoms = 5;
    h_amino_acids[AA_ALA].n_sidechain_atoms = 2;
    h_amino_acids[AA_ALA].atom_order[0] = ATOM_N;
    h_amino_acids[AA_ALA].atom_order[1] = ATOM_CA;
    h_amino_acids[AA_ALA].atom_order[2] = ATOM_C;
    h_amino_acids[AA_ALA].atom_order[3] = ATOM_O;
    h_amino_acids[AA_ALA].atom_order[4] = ATOM_CB;
    h_amino_acids[AA_ALA].alt_atom_order[0] = ATOM_N;
    h_amino_acids[AA_ALA].alt_atom_order[1] = ATOM_CA;
    h_amino_acids[AA_ALA].alt_atom_order[2] = ATOM_C;
    h_amino_acids[AA_ALA].alt_atom_order[3] = ATOM_CB;
    h_amino_acids[AA_ALA].alt_atom_order[4] = ATOM_O;
    SET_SC_ATOM(AA_ALA, 0, ATOM_O, ATOM_N, ATOM_CA, ATOM_C, 1.23f, 120.31f);
    SET_SC_ATOM(AA_ALA, 1, ATOM_CB, ATOM_O, ATOM_C, ATOM_CA, 1.52f, 110.852f);

    // ARG - Arginine
    h_amino_acids[AA_ARG].n_atoms = 11;
    h_amino_acids[AA_ARG].n_sidechain_atoms = 8;
    h_amino_acids[AA_ARG].atom_order[0] = ATOM_N;
    h_amino_acids[AA_ARG].atom_order[1] = ATOM_CA;
    h_amino_acids[AA_ARG].atom_order[2] = ATOM_C;
    h_amino_acids[AA_ARG].atom_order[3] = ATOM_O;
    h_amino_acids[AA_ARG].atom_order[4] = ATOM_CB;
    h_amino_acids[AA_ARG].atom_order[5] = ATOM_CG;
    h_amino_acids[AA_ARG].atom_order[6] = ATOM_CD;
    h_amino_acids[AA_ARG].atom_order[7] = ATOM_NE;
    h_amino_acids[AA_ARG].atom_order[8] = ATOM_CZ;
    h_amino_acids[AA_ARG].atom_order[9] = ATOM_NH1;
    h_amino_acids[AA_ARG].atom_order[10] = ATOM_NH2;
    h_amino_acids[AA_ARG].alt_atom_order[0] = ATOM_N;
    h_amino_acids[AA_ARG].alt_atom_order[1] = ATOM_CA;
    h_amino_acids[AA_ARG].alt_atom_order[2] = ATOM_C;
    h_amino_acids[AA_ARG].alt_atom_order[3] = ATOM_CB;
    h_amino_acids[AA_ARG].alt_atom_order[4] = ATOM_O;
    h_amino_acids[AA_ARG].alt_atom_order[5] = ATOM_CG;
    h_amino_acids[AA_ARG].alt_atom_order[6] = ATOM_CD;
    h_amino_acids[AA_ARG].alt_atom_order[7] = ATOM_NE;
    h_amino_acids[AA_ARG].alt_atom_order[8] = ATOM_NH1;
    h_amino_acids[AA_ARG].alt_atom_order[9] = ATOM_NH2;
    h_amino_acids[AA_ARG].alt_atom_order[10] = ATOM_CZ;
    SET_SC_ATOM(AA_ARG, 0, ATOM_O, ATOM_N, ATOM_CA, ATOM_C, 1.23f, 119.745f);
    SET_SC_ATOM(AA_ARG, 1, ATOM_CB, ATOM_O, ATOM_C, ATOM_CA, 1.53f, 110.579f);
    SET_SC_ATOM(AA_ARG, 2, ATOM_CG, ATOM_N, ATOM_CA, ATOM_CB, 1.53f, 113.233f);
    SET_SC_ATOM(AA_ARG, 3, ATOM_CD, ATOM_CA, ATOM_CB, ATOM_CG, 1.52f, 110.787f);
    SET_SC_ATOM(AA_ARG, 4, ATOM_NE, ATOM_CB, ATOM_CG, ATOM_CD, 1.46f, 111.919f);
    SET_SC_ATOM(AA_ARG, 5, ATOM_CZ, ATOM_CG, ATOM_CD, ATOM_NE, 1.32f, 125.192f);
    SET_SC_ATOM(AA_ARG, 6, ATOM_NH1, ATOM_CD, ATOM_NE, ATOM_CZ, 1.31f, 120.077f);
    SET_SC_ATOM(AA_ARG, 7, ATOM_NH2, ATOM_CD, ATOM_NE, ATOM_CZ, 1.31f, 120.077f);

    // ASN - Asparagine
    h_amino_acids[AA_ASN].n_atoms = 8;
    h_amino_acids[AA_ASN].n_sidechain_atoms = 5;
    h_amino_acids[AA_ASN].atom_order[0] = ATOM_N;
    h_amino_acids[AA_ASN].atom_order[1] = ATOM_CA;
    h_amino_acids[AA_ASN].atom_order[2] = ATOM_C;
    h_amino_acids[AA_ASN].atom_order[3] = ATOM_O;
    h_amino_acids[AA_ASN].atom_order[4] = ATOM_CB;
    h_amino_acids[AA_ASN].atom_order[5] = ATOM_CG;
    h_amino_acids[AA_ASN].atom_order[6] = ATOM_OD1;
    h_amino_acids[AA_ASN].atom_order[7] = ATOM_ND2;
    h_amino_acids[AA_ASN].alt_atom_order[0] = ATOM_N;
    h_amino_acids[AA_ASN].alt_atom_order[1] = ATOM_CA;
    h_amino_acids[AA_ASN].alt_atom_order[2] = ATOM_C;
    h_amino_acids[AA_ASN].alt_atom_order[3] = ATOM_CB;
    h_amino_acids[AA_ASN].alt_atom_order[4] = ATOM_O;
    h_amino_acids[AA_ASN].alt_atom_order[5] = ATOM_CG;
    h_amino_acids[AA_ASN].alt_atom_order[6] = ATOM_ND2;
    h_amino_acids[AA_ASN].alt_atom_order[7] = ATOM_OD1;
    SET_SC_ATOM(AA_ASN, 0, ATOM_O, ATOM_N, ATOM_CA, ATOM_C, 1.23f, 120.313f);
    SET_SC_ATOM(AA_ASN, 1, ATOM_CB, ATOM_O, ATOM_C, ATOM_CA, 1.52f, 110.852f);
    SET_SC_ATOM(AA_ASN, 2, ATOM_CG, ATOM_N, ATOM_CA, ATOM_CB, 1.52f, 113.232f);
    SET_SC_ATOM(AA_ASN, 3, ATOM_OD1, ATOM_CA, ATOM_CB, ATOM_CG, 1.23f, 120.85f);
    SET_SC_ATOM(AA_ASN, 4, ATOM_ND2, ATOM_CA, ATOM_CB, ATOM_CG, 1.325f, 116.48f);

    // ASP - Aspartic acid
    h_amino_acids[AA_ASP].n_atoms = 8;
    h_amino_acids[AA_ASP].n_sidechain_atoms = 5;
    h_amino_acids[AA_ASP].atom_order[0] = ATOM_N;
    h_amino_acids[AA_ASP].atom_order[1] = ATOM_CA;
    h_amino_acids[AA_ASP].atom_order[2] = ATOM_C;
    h_amino_acids[AA_ASP].atom_order[3] = ATOM_O;
    h_amino_acids[AA_ASP].atom_order[4] = ATOM_CB;
    h_amino_acids[AA_ASP].atom_order[5] = ATOM_CG;
    h_amino_acids[AA_ASP].atom_order[6] = ATOM_OD1;
    h_amino_acids[AA_ASP].atom_order[7] = ATOM_OD2;
    h_amino_acids[AA_ASP].alt_atom_order[0] = ATOM_N;
    h_amino_acids[AA_ASP].alt_atom_order[1] = ATOM_CA;
    h_amino_acids[AA_ASP].alt_atom_order[2] = ATOM_C;
    h_amino_acids[AA_ASP].alt_atom_order[3] = ATOM_CB;
    h_amino_acids[AA_ASP].alt_atom_order[4] = ATOM_O;
    h_amino_acids[AA_ASP].alt_atom_order[5] = ATOM_CG;
    h_amino_acids[AA_ASP].alt_atom_order[6] = ATOM_OD1;
    h_amino_acids[AA_ASP].alt_atom_order[7] = ATOM_OD2;
    SET_SC_ATOM(AA_ASP, 0, ATOM_O, ATOM_N, ATOM_CA, ATOM_C, 1.23f, 121.051f);
    SET_SC_ATOM(AA_ASP, 1, ATOM_CB, ATOM_O, ATOM_C, ATOM_CA, 1.53f, 110.871f);
    SET_SC_ATOM(AA_ASP, 2, ATOM_CG, ATOM_N, ATOM_CA, ATOM_CB, 1.52f, 113.232f);
    SET_SC_ATOM(AA_ASP, 3, ATOM_OD1, ATOM_CA, ATOM_CB, ATOM_CG, 1.248f, 118.344f);
    SET_SC_ATOM(AA_ASP, 4, ATOM_OD2, ATOM_CA, ATOM_CB, ATOM_CG, 1.248f, 118.344f);

    // CYS - Cysteine
    h_amino_acids[AA_CYS].n_atoms = 6;
    h_amino_acids[AA_CYS].n_sidechain_atoms = 3;
    h_amino_acids[AA_CYS].atom_order[0] = ATOM_N;
    h_amino_acids[AA_CYS].atom_order[1] = ATOM_CA;
    h_amino_acids[AA_CYS].atom_order[2] = ATOM_C;
    h_amino_acids[AA_CYS].atom_order[3] = ATOM_O;
    h_amino_acids[AA_CYS].atom_order[4] = ATOM_CB;
    h_amino_acids[AA_CYS].atom_order[5] = ATOM_SG;
    h_amino_acids[AA_CYS].alt_atom_order[0] = ATOM_N;
    h_amino_acids[AA_CYS].alt_atom_order[1] = ATOM_CA;
    h_amino_acids[AA_CYS].alt_atom_order[2] = ATOM_C;
    h_amino_acids[AA_CYS].alt_atom_order[3] = ATOM_CB;
    h_amino_acids[AA_CYS].alt_atom_order[4] = ATOM_O;
    h_amino_acids[AA_CYS].alt_atom_order[5] = ATOM_SG;
    SET_SC_ATOM(AA_CYS, 0, ATOM_O, ATOM_N, ATOM_CA, ATOM_C, 1.23f, 120.063f);
    SET_SC_ATOM(AA_CYS, 1, ATOM_CB, ATOM_O, ATOM_C, ATOM_CA, 1.53f, 111.078f);
    SET_SC_ATOM(AA_CYS, 2, ATOM_SG, ATOM_N, ATOM_CA, ATOM_CB, 1.8f, 113.817f);

    // GLN - Glutamine
    h_amino_acids[AA_GLN].n_atoms = 9;
    h_amino_acids[AA_GLN].n_sidechain_atoms = 6;
    h_amino_acids[AA_GLN].atom_order[0] = ATOM_N;
    h_amino_acids[AA_GLN].atom_order[1] = ATOM_CA;
    h_amino_acids[AA_GLN].atom_order[2] = ATOM_C;
    h_amino_acids[AA_GLN].atom_order[3] = ATOM_O;
    h_amino_acids[AA_GLN].atom_order[4] = ATOM_CB;
    h_amino_acids[AA_GLN].atom_order[5] = ATOM_CG;
    h_amino_acids[AA_GLN].atom_order[6] = ATOM_CD;
    h_amino_acids[AA_GLN].atom_order[7] = ATOM_OE1;
    h_amino_acids[AA_GLN].atom_order[8] = ATOM_NE2;
    h_amino_acids[AA_GLN].alt_atom_order[0] = ATOM_N;
    h_amino_acids[AA_GLN].alt_atom_order[1] = ATOM_CA;
    h_amino_acids[AA_GLN].alt_atom_order[2] = ATOM_C;
    h_amino_acids[AA_GLN].alt_atom_order[3] = ATOM_CB;
    h_amino_acids[AA_GLN].alt_atom_order[4] = ATOM_O;
    h_amino_acids[AA_GLN].alt_atom_order[5] = ATOM_CG;
    h_amino_acids[AA_GLN].alt_atom_order[6] = ATOM_CD;
    h_amino_acids[AA_GLN].alt_atom_order[7] = ATOM_NE2;
    h_amino_acids[AA_GLN].alt_atom_order[8] = ATOM_OE1;
    SET_SC_ATOM(AA_GLN, 0, ATOM_O, ATOM_N, ATOM_CA, ATOM_C, 1.23f, 120.211f);
    SET_SC_ATOM(AA_GLN, 1, ATOM_CB, ATOM_O, ATOM_C, ATOM_CA, 1.53f, 109.5f);
    SET_SC_ATOM(AA_GLN, 2, ATOM_CG, ATOM_N, ATOM_CA, ATOM_CB, 1.52f, 113.292f);
    SET_SC_ATOM(AA_GLN, 3, ATOM_CD, ATOM_CA, ATOM_CB, ATOM_CG, 1.52f, 112.811f);
    SET_SC_ATOM(AA_GLN, 4, ATOM_OE1, ATOM_CB, ATOM_CG, ATOM_CD, 1.23f, 121.844f);
    SET_SC_ATOM(AA_GLN, 5, ATOM_NE2, ATOM_CB, ATOM_CG, ATOM_CD, 1.32f, 116.50f);

    // GLU - Glutamic acid
    h_amino_acids[AA_GLU].n_atoms = 9;
    h_amino_acids[AA_GLU].n_sidechain_atoms = 6;
    h_amino_acids[AA_GLU].atom_order[0] = ATOM_N;
    h_amino_acids[AA_GLU].atom_order[1] = ATOM_CA;
    h_amino_acids[AA_GLU].atom_order[2] = ATOM_C;
    h_amino_acids[AA_GLU].atom_order[3] = ATOM_O;
    h_amino_acids[AA_GLU].atom_order[4] = ATOM_CB;
    h_amino_acids[AA_GLU].atom_order[5] = ATOM_CG;
    h_amino_acids[AA_GLU].atom_order[6] = ATOM_CD;
    h_amino_acids[AA_GLU].atom_order[7] = ATOM_OE1;
    h_amino_acids[AA_GLU].atom_order[8] = ATOM_OE2;
    h_amino_acids[AA_GLU].alt_atom_order[0] = ATOM_N;
    h_amino_acids[AA_GLU].alt_atom_order[1] = ATOM_CA;
    h_amino_acids[AA_GLU].alt_atom_order[2] = ATOM_C;
    h_amino_acids[AA_GLU].alt_atom_order[3] = ATOM_CB;
    h_amino_acids[AA_GLU].alt_atom_order[4] = ATOM_O;
    h_amino_acids[AA_GLU].alt_atom_order[5] = ATOM_CG;
    h_amino_acids[AA_GLU].alt_atom_order[6] = ATOM_CD;
    h_amino_acids[AA_GLU].alt_atom_order[7] = ATOM_OE1;
    h_amino_acids[AA_GLU].alt_atom_order[8] = ATOM_OE2;
    SET_SC_ATOM(AA_GLU, 0, ATOM_O, ATOM_N, ATOM_CA, ATOM_C, 1.23f, 120.594f);
    SET_SC_ATOM(AA_GLU, 1, ATOM_CB, ATOM_O, ATOM_C, ATOM_CA, 1.53f, 110.538f);
    SET_SC_ATOM(AA_GLU, 2, ATOM_CG, ATOM_N, ATOM_CA, ATOM_CB, 1.52f, 113.82f);
    SET_SC_ATOM(AA_GLU, 3, ATOM_CD, ATOM_CA, ATOM_CB, ATOM_CG, 1.52f, 112.912f);
    SET_SC_ATOM(AA_GLU, 4, ATOM_OE1, ATOM_CB, ATOM_CG, ATOM_CD, 1.25f, 118.479f);
    SET_SC_ATOM(AA_GLU, 5, ATOM_OE2, ATOM_CB, ATOM_CG, ATOM_CD, 1.25f, 118.479f);

    // GLY - Glycine
    h_amino_acids[AA_GLY].n_atoms = 4;
    h_amino_acids[AA_GLY].n_sidechain_atoms = 1;
    h_amino_acids[AA_GLY].atom_order[0] = ATOM_N;
    h_amino_acids[AA_GLY].atom_order[1] = ATOM_CA;
    h_amino_acids[AA_GLY].atom_order[2] = ATOM_C;
    h_amino_acids[AA_GLY].atom_order[3] = ATOM_O;
    h_amino_acids[AA_GLY].alt_atom_order[0] = ATOM_N;
    h_amino_acids[AA_GLY].alt_atom_order[1] = ATOM_CA;
    h_amino_acids[AA_GLY].alt_atom_order[2] = ATOM_C;
    h_amino_acids[AA_GLY].alt_atom_order[3] = ATOM_O;
    SET_SC_ATOM(AA_GLY, 0, ATOM_O, ATOM_N, ATOM_CA, ATOM_C, 1.23f, 120.522f);

    // HIS - Histidine
    h_amino_acids[AA_HIS].n_atoms = 10;
    h_amino_acids[AA_HIS].n_sidechain_atoms = 7;
    h_amino_acids[AA_HIS].atom_order[0] = ATOM_N;
    h_amino_acids[AA_HIS].atom_order[1] = ATOM_CA;
    h_amino_acids[AA_HIS].atom_order[2] = ATOM_C;
    h_amino_acids[AA_HIS].atom_order[3] = ATOM_O;
    h_amino_acids[AA_HIS].atom_order[4] = ATOM_CB;
    h_amino_acids[AA_HIS].atom_order[5] = ATOM_CG;
    h_amino_acids[AA_HIS].atom_order[6] = ATOM_ND1;
    h_amino_acids[AA_HIS].atom_order[7] = ATOM_CD2;
    h_amino_acids[AA_HIS].atom_order[8] = ATOM_CE1;
    h_amino_acids[AA_HIS].atom_order[9] = ATOM_NE2;
    h_amino_acids[AA_HIS].alt_atom_order[0] = ATOM_N;
    h_amino_acids[AA_HIS].alt_atom_order[1] = ATOM_CA;
    h_amino_acids[AA_HIS].alt_atom_order[2] = ATOM_C;
    h_amino_acids[AA_HIS].alt_atom_order[3] = ATOM_CB;
    h_amino_acids[AA_HIS].alt_atom_order[4] = ATOM_O;
    h_amino_acids[AA_HIS].alt_atom_order[5] = ATOM_CG;
    h_amino_acids[AA_HIS].alt_atom_order[6] = ATOM_CD2;
    h_amino_acids[AA_HIS].alt_atom_order[7] = ATOM_ND1;
    h_amino_acids[AA_HIS].alt_atom_order[8] = ATOM_CE1;
    h_amino_acids[AA_HIS].alt_atom_order[9] = ATOM_NE2;
    SET_SC_ATOM(AA_HIS, 0, ATOM_O, ATOM_N, ATOM_CA, ATOM_C, 1.23f, 120.548f);
    SET_SC_ATOM(AA_HIS, 1, ATOM_CB, ATOM_O, ATOM_C, ATOM_CA, 1.53f, 111.329f);
    SET_SC_ATOM(AA_HIS, 2, ATOM_CG, ATOM_N, ATOM_CA, ATOM_CB, 1.5f, 113.468f);
    SET_SC_ATOM(AA_HIS, 3, ATOM_ND1, ATOM_CA, ATOM_CB, ATOM_CG, 1.38f, 122.85f);
    SET_SC_ATOM(AA_HIS, 4, ATOM_CD2, ATOM_CA, ATOM_CB, ATOM_CG, 1.36f, 130.61f);
    SET_SC_ATOM(AA_HIS, 5, ATOM_CE1, ATOM_CB, ATOM_CG, ATOM_ND1, 1.33f, 108.589f);
    SET_SC_ATOM(AA_HIS, 6, ATOM_NE2, ATOM_CB, ATOM_CG, ATOM_CD2, 1.38f, 107.439f);

    // ILE - Isoleucine
    h_amino_acids[AA_ILE].n_atoms = 8;
    h_amino_acids[AA_ILE].n_sidechain_atoms = 5;
    h_amino_acids[AA_ILE].atom_order[0] = ATOM_N;
    h_amino_acids[AA_ILE].atom_order[1] = ATOM_CA;
    h_amino_acids[AA_ILE].atom_order[2] = ATOM_C;
    h_amino_acids[AA_ILE].atom_order[3] = ATOM_O;
    h_amino_acids[AA_ILE].atom_order[4] = ATOM_CB;
    h_amino_acids[AA_ILE].atom_order[5] = ATOM_CG1;
    h_amino_acids[AA_ILE].atom_order[6] = ATOM_CG2;
    h_amino_acids[AA_ILE].atom_order[7] = ATOM_CD1;
    h_amino_acids[AA_ILE].alt_atom_order[0] = ATOM_N;
    h_amino_acids[AA_ILE].alt_atom_order[1] = ATOM_CA;
    h_amino_acids[AA_ILE].alt_atom_order[2] = ATOM_C;
    h_amino_acids[AA_ILE].alt_atom_order[3] = ATOM_CB;
    h_amino_acids[AA_ILE].alt_atom_order[4] = ATOM_O;
    h_amino_acids[AA_ILE].alt_atom_order[5] = ATOM_CG1;
    h_amino_acids[AA_ILE].alt_atom_order[6] = ATOM_CG2;
    h_amino_acids[AA_ILE].alt_atom_order[7] = ATOM_CD1;
    SET_SC_ATOM(AA_ILE, 0, ATOM_O, ATOM_N, ATOM_CA, ATOM_C, 1.235f, 120.393f);
    SET_SC_ATOM(AA_ILE, 1, ATOM_CB, ATOM_O, ATOM_C, ATOM_CA, 1.54f, 111.983f);
    SET_SC_ATOM(AA_ILE, 2, ATOM_CG1, ATOM_N, ATOM_CA, ATOM_CB, 1.53f, 110.5f);
    SET_SC_ATOM(AA_ILE, 3, ATOM_CG2, ATOM_N, ATOM_CA, ATOM_CB, 1.52f, 110.5f);
    SET_SC_ATOM(AA_ILE, 4, ATOM_CD1, ATOM_CA, ATOM_CB, ATOM_CG1, 1.51f, 113.97f);

    // LEU - Leucine
    h_amino_acids[AA_LEU].n_atoms = 8;
    h_amino_acids[AA_LEU].n_sidechain_atoms = 5;
    h_amino_acids[AA_LEU].atom_order[0] = ATOM_N;
    h_amino_acids[AA_LEU].atom_order[1] = ATOM_CA;
    h_amino_acids[AA_LEU].atom_order[2] = ATOM_C;
    h_amino_acids[AA_LEU].atom_order[3] = ATOM_O;
    h_amino_acids[AA_LEU].atom_order[4] = ATOM_CB;
    h_amino_acids[AA_LEU].atom_order[5] = ATOM_CG;
    h_amino_acids[AA_LEU].atom_order[6] = ATOM_CD1;
    h_amino_acids[AA_LEU].atom_order[7] = ATOM_CD2;
    h_amino_acids[AA_LEU].alt_atom_order[0] = ATOM_N;
    h_amino_acids[AA_LEU].alt_atom_order[1] = ATOM_CA;
    h_amino_acids[AA_LEU].alt_atom_order[2] = ATOM_C;
    h_amino_acids[AA_LEU].alt_atom_order[3] = ATOM_CB;
    h_amino_acids[AA_LEU].alt_atom_order[4] = ATOM_O;
    h_amino_acids[AA_LEU].alt_atom_order[5] = ATOM_CG;
    h_amino_acids[AA_LEU].alt_atom_order[6] = ATOM_CD1;
    h_amino_acids[AA_LEU].alt_atom_order[7] = ATOM_CD2;
    SET_SC_ATOM(AA_LEU, 0, ATOM_O, ATOM_N, ATOM_CA, ATOM_C, 1.235f, 120.211f);
    SET_SC_ATOM(AA_LEU, 1, ATOM_CB, ATOM_O, ATOM_C, ATOM_CA, 1.53f, 110.418f);
    SET_SC_ATOM(AA_LEU, 2, ATOM_CG, ATOM_N, ATOM_CA, ATOM_CB, 1.53f, 116.10f);
    SET_SC_ATOM(AA_LEU, 3, ATOM_CD1, ATOM_CA, ATOM_CB, ATOM_CG, 1.52f, 110.58f);
    SET_SC_ATOM(AA_LEU, 4, ATOM_CD2, ATOM_CA, ATOM_CB, ATOM_CG, 1.52f, 110.58f);

    // LYS - Lysine
    h_amino_acids[AA_LYS].n_atoms = 9;
    h_amino_acids[AA_LYS].n_sidechain_atoms = 6;
    h_amino_acids[AA_LYS].atom_order[0] = ATOM_N;
    h_amino_acids[AA_LYS].atom_order[1] = ATOM_CA;
    h_amino_acids[AA_LYS].atom_order[2] = ATOM_C;
    h_amino_acids[AA_LYS].atom_order[3] = ATOM_O;
    h_amino_acids[AA_LYS].atom_order[4] = ATOM_CB;
    h_amino_acids[AA_LYS].atom_order[5] = ATOM_CG;
    h_amino_acids[AA_LYS].atom_order[6] = ATOM_CD;
    h_amino_acids[AA_LYS].atom_order[7] = ATOM_CE;
    h_amino_acids[AA_LYS].atom_order[8] = ATOM_NZ;
    h_amino_acids[AA_LYS].alt_atom_order[0] = ATOM_N;
    h_amino_acids[AA_LYS].alt_atom_order[1] = ATOM_CA;
    h_amino_acids[AA_LYS].alt_atom_order[2] = ATOM_C;
    h_amino_acids[AA_LYS].alt_atom_order[3] = ATOM_CB;
    h_amino_acids[AA_LYS].alt_atom_order[4] = ATOM_O;
    h_amino_acids[AA_LYS].alt_atom_order[5] = ATOM_CG;
    h_amino_acids[AA_LYS].alt_atom_order[6] = ATOM_CD;
    h_amino_acids[AA_LYS].alt_atom_order[7] = ATOM_CE;
    h_amino_acids[AA_LYS].alt_atom_order[8] = ATOM_NZ;
    SET_SC_ATOM(AA_LYS, 0, ATOM_O, ATOM_N, ATOM_CA, ATOM_C, 1.23f, 120.54f);
    SET_SC_ATOM(AA_LYS, 1, ATOM_CB, ATOM_O, ATOM_C, ATOM_CA, 1.53f, 109.5f);
    SET_SC_ATOM(AA_LYS, 2, ATOM_CG, ATOM_N, ATOM_CA, ATOM_CB, 1.52f, 113.83f);
    SET_SC_ATOM(AA_LYS, 3, ATOM_CD, ATOM_CA, ATOM_CB, ATOM_CG, 1.52f, 111.79f);
    SET_SC_ATOM(AA_LYS, 4, ATOM_CE, ATOM_CB, ATOM_CG, ATOM_CD, 1.52f, 111.79f);
    SET_SC_ATOM(AA_LYS, 5, ATOM_NZ, ATOM_CG, ATOM_CD, ATOM_CE, 1.49f, 112.25f);

    // MET - Methionine
    h_amino_acids[AA_MET].n_atoms = 8;
    h_amino_acids[AA_MET].n_sidechain_atoms = 5;
    h_amino_acids[AA_MET].atom_order[0] = ATOM_N;
    h_amino_acids[AA_MET].atom_order[1] = ATOM_CA;
    h_amino_acids[AA_MET].atom_order[2] = ATOM_C;
    h_amino_acids[AA_MET].atom_order[3] = ATOM_O;
    h_amino_acids[AA_MET].atom_order[4] = ATOM_CB;
    h_amino_acids[AA_MET].atom_order[5] = ATOM_CG;
    h_amino_acids[AA_MET].atom_order[6] = ATOM_SD;
    h_amino_acids[AA_MET].atom_order[7] = ATOM_CE;
    h_amino_acids[AA_MET].alt_atom_order[0] = ATOM_N;
    h_amino_acids[AA_MET].alt_atom_order[1] = ATOM_CA;
    h_amino_acids[AA_MET].alt_atom_order[2] = ATOM_C;
    h_amino_acids[AA_MET].alt_atom_order[3] = ATOM_CB;
    h_amino_acids[AA_MET].alt_atom_order[4] = ATOM_O;
    h_amino_acids[AA_MET].alt_atom_order[5] = ATOM_CG;
    h_amino_acids[AA_MET].alt_atom_order[6] = ATOM_SD;
    h_amino_acids[AA_MET].alt_atom_order[7] = ATOM_CE;
    SET_SC_ATOM(AA_MET, 0, ATOM_O, ATOM_N, ATOM_CA, ATOM_C, 1.23f, 120.148f);
    SET_SC_ATOM(AA_MET, 1, ATOM_CB, ATOM_O, ATOM_C, ATOM_CA, 1.53f, 110.833f);
    SET_SC_ATOM(AA_MET, 2, ATOM_CG, ATOM_N, ATOM_CA, ATOM_CB, 1.52f, 113.68f);
    SET_SC_ATOM(AA_MET, 3, ATOM_SD, ATOM_CA, ATOM_CB, ATOM_CG, 1.8f, 112.773f);
    SET_SC_ATOM(AA_MET, 4, ATOM_CE, ATOM_CB, ATOM_CG, ATOM_SD, 1.79f, 100.61f);

    // PHE - Phenylalanine
    h_amino_acids[AA_PHE].n_atoms = 11;
    h_amino_acids[AA_PHE].n_sidechain_atoms = 8;
    h_amino_acids[AA_PHE].atom_order[0] = ATOM_N;
    h_amino_acids[AA_PHE].atom_order[1] = ATOM_CA;
    h_amino_acids[AA_PHE].atom_order[2] = ATOM_C;
    h_amino_acids[AA_PHE].atom_order[3] = ATOM_O;
    h_amino_acids[AA_PHE].atom_order[4] = ATOM_CB;
    h_amino_acids[AA_PHE].atom_order[5] = ATOM_CG;
    h_amino_acids[AA_PHE].atom_order[6] = ATOM_CD1;
    h_amino_acids[AA_PHE].atom_order[7] = ATOM_CD2;
    h_amino_acids[AA_PHE].atom_order[8] = ATOM_CE1;
    h_amino_acids[AA_PHE].atom_order[9] = ATOM_CE2;
    h_amino_acids[AA_PHE].atom_order[10] = ATOM_CZ;
    h_amino_acids[AA_PHE].alt_atom_order[0] = ATOM_N;
    h_amino_acids[AA_PHE].alt_atom_order[1] = ATOM_CA;
    h_amino_acids[AA_PHE].alt_atom_order[2] = ATOM_C;
    h_amino_acids[AA_PHE].alt_atom_order[3] = ATOM_CB;
    h_amino_acids[AA_PHE].alt_atom_order[4] = ATOM_O;
    h_amino_acids[AA_PHE].alt_atom_order[5] = ATOM_CG;
    h_amino_acids[AA_PHE].alt_atom_order[6] = ATOM_CD1;
    h_amino_acids[AA_PHE].alt_atom_order[7] = ATOM_CD2;
    h_amino_acids[AA_PHE].alt_atom_order[8] = ATOM_CE1;
    h_amino_acids[AA_PHE].alt_atom_order[9] = ATOM_CE2;
    h_amino_acids[AA_PHE].alt_atom_order[10] = ATOM_CZ;
    SET_SC_ATOM(AA_PHE, 0, ATOM_O, ATOM_N, ATOM_CA, ATOM_C, 1.23f, 120.283f);
    SET_SC_ATOM(AA_PHE, 1, ATOM_CB, ATOM_O, ATOM_C, ATOM_CA, 1.53f, 110.846f);
    SET_SC_ATOM(AA_PHE, 2, ATOM_CG, ATOM_N, ATOM_CA, ATOM_CB, 1.51f, 114.0f);
    SET_SC_ATOM(AA_PHE, 3, ATOM_CD1, ATOM_CA, ATOM_CB, ATOM_CG, 1.385f, 120.0f);
    SET_SC_ATOM(AA_PHE, 4, ATOM_CD2, ATOM_CA, ATOM_CB, ATOM_CG, 1.385f, 120.0f);
    SET_SC_ATOM(AA_PHE, 5, ATOM_CE1, ATOM_CB, ATOM_CG, ATOM_CD1, 1.385f, 120.0f);
    SET_SC_ATOM(AA_PHE, 6, ATOM_CE2, ATOM_CB, ATOM_CG, ATOM_CD2, 1.385f, 120.0f);
    SET_SC_ATOM(AA_PHE, 7, ATOM_CZ, ATOM_CG, ATOM_CD1, ATOM_CE1, 1.385f, 120.0f);

    // PRO - Proline
    h_amino_acids[AA_PRO].n_atoms = 7;
    h_amino_acids[AA_PRO].n_sidechain_atoms = 4;
    h_amino_acids[AA_PRO].atom_order[0] = ATOM_N;
    h_amino_acids[AA_PRO].atom_order[1] = ATOM_CA;
    h_amino_acids[AA_PRO].atom_order[2] = ATOM_C;
    h_amino_acids[AA_PRO].atom_order[3] = ATOM_O;
    h_amino_acids[AA_PRO].atom_order[4] = ATOM_CB;
    h_amino_acids[AA_PRO].atom_order[5] = ATOM_CG;
    h_amino_acids[AA_PRO].atom_order[6] = ATOM_CD;
    h_amino_acids[AA_PRO].alt_atom_order[0] = ATOM_N;
    h_amino_acids[AA_PRO].alt_atom_order[1] = ATOM_CA;
    h_amino_acids[AA_PRO].alt_atom_order[2] = ATOM_C;
    h_amino_acids[AA_PRO].alt_atom_order[3] = ATOM_CB;
    h_amino_acids[AA_PRO].alt_atom_order[4] = ATOM_O;
    h_amino_acids[AA_PRO].alt_atom_order[5] = ATOM_CG;
    h_amino_acids[AA_PRO].alt_atom_order[6] = ATOM_CD;
    SET_SC_ATOM(AA_PRO, 0, ATOM_O, ATOM_N, ATOM_CA, ATOM_C, 1.23f, 120.6f);
    SET_SC_ATOM(AA_PRO, 1, ATOM_CB, ATOM_O, ATOM_C, ATOM_CA, 1.53f, 111.372f);
    SET_SC_ATOM(AA_PRO, 2, ATOM_CG, ATOM_N, ATOM_CA, ATOM_CB, 1.49f, 104.21f);
    SET_SC_ATOM(AA_PRO, 3, ATOM_CD, ATOM_CA, ATOM_CB, ATOM_CG, 1.50f, 105.0f);

    // SER - Serine
    h_amino_acids[AA_SER].n_atoms = 6;
    h_amino_acids[AA_SER].n_sidechain_atoms = 3;
    h_amino_acids[AA_SER].atom_order[0] = ATOM_N;
    h_amino_acids[AA_SER].atom_order[1] = ATOM_CA;
    h_amino_acids[AA_SER].atom_order[2] = ATOM_C;
    h_amino_acids[AA_SER].atom_order[3] = ATOM_O;
    h_amino_acids[AA_SER].atom_order[4] = ATOM_CB;
    h_amino_acids[AA_SER].atom_order[5] = ATOM_OG;
    h_amino_acids[AA_SER].alt_atom_order[0] = ATOM_N;
    h_amino_acids[AA_SER].alt_atom_order[1] = ATOM_CA;
    h_amino_acids[AA_SER].alt_atom_order[2] = ATOM_C;
    h_amino_acids[AA_SER].alt_atom_order[3] = ATOM_CB;
    h_amino_acids[AA_SER].alt_atom_order[4] = ATOM_O;
    h_amino_acids[AA_SER].alt_atom_order[5] = ATOM_OG;
    SET_SC_ATOM(AA_SER, 0, ATOM_O, ATOM_N, ATOM_CA, ATOM_C, 1.23f, 120.475f);
    SET_SC_ATOM(AA_SER, 1, ATOM_CB, ATOM_O, ATOM_C, ATOM_CA, 1.53f, 110.248f);
    SET_SC_ATOM(AA_SER, 2, ATOM_OG, ATOM_N, ATOM_CA, ATOM_CB, 1.417f, 111.132f);

    // THR - Threonine
    h_amino_acids[AA_THR].n_atoms = 7;
    h_amino_acids[AA_THR].n_sidechain_atoms = 4;
    h_amino_acids[AA_THR].atom_order[0] = ATOM_N;
    h_amino_acids[AA_THR].atom_order[1] = ATOM_CA;
    h_amino_acids[AA_THR].atom_order[2] = ATOM_C;
    h_amino_acids[AA_THR].atom_order[3] = ATOM_O;
    h_amino_acids[AA_THR].atom_order[4] = ATOM_CB;
    h_amino_acids[AA_THR].atom_order[5] = ATOM_OG1;
    h_amino_acids[AA_THR].atom_order[6] = ATOM_CG2;
    h_amino_acids[AA_THR].alt_atom_order[0] = ATOM_N;
    h_amino_acids[AA_THR].alt_atom_order[1] = ATOM_CA;
    h_amino_acids[AA_THR].alt_atom_order[2] = ATOM_C;
    h_amino_acids[AA_THR].alt_atom_order[3] = ATOM_CB;
    h_amino_acids[AA_THR].alt_atom_order[4] = ATOM_O;
    h_amino_acids[AA_THR].alt_atom_order[5] = ATOM_CG2;
    h_amino_acids[AA_THR].alt_atom_order[6] = ATOM_OG1;
    SET_SC_ATOM(AA_THR, 0, ATOM_O, ATOM_N, ATOM_CA, ATOM_C, 1.23f, 120.252f);
    SET_SC_ATOM(AA_THR, 1, ATOM_CB, ATOM_O, ATOM_C, ATOM_CA, 1.53f, 110.075f);
    SET_SC_ATOM(AA_THR, 2, ATOM_OG1, ATOM_N, ATOM_CA, ATOM_CB, 1.43f, 109.442f);
    SET_SC_ATOM(AA_THR, 3, ATOM_CG2, ATOM_N, ATOM_CA, ATOM_CB, 1.52f, 111.457f);

    // TRP - Tryptophan
    h_amino_acids[AA_TRP].n_atoms = 14;
    h_amino_acids[AA_TRP].n_sidechain_atoms = 11;
    h_amino_acids[AA_TRP].atom_order[0] = ATOM_N;
    h_amino_acids[AA_TRP].atom_order[1] = ATOM_CA;
    h_amino_acids[AA_TRP].atom_order[2] = ATOM_C;
    h_amino_acids[AA_TRP].atom_order[3] = ATOM_O;
    h_amino_acids[AA_TRP].atom_order[4] = ATOM_CB;
    h_amino_acids[AA_TRP].atom_order[5] = ATOM_CG;
    h_amino_acids[AA_TRP].atom_order[6] = ATOM_CD1;
    h_amino_acids[AA_TRP].atom_order[7] = ATOM_CD2;
    h_amino_acids[AA_TRP].atom_order[8] = ATOM_NE1;
    h_amino_acids[AA_TRP].atom_order[9] = ATOM_CE2;
    h_amino_acids[AA_TRP].atom_order[10] = ATOM_CE3;
    h_amino_acids[AA_TRP].atom_order[11] = ATOM_CZ2;
    h_amino_acids[AA_TRP].atom_order[12] = ATOM_CZ3;
    h_amino_acids[AA_TRP].atom_order[13] = ATOM_CH2;
    h_amino_acids[AA_TRP].alt_atom_order[0] = ATOM_N;
    h_amino_acids[AA_TRP].alt_atom_order[1] = ATOM_CA;
    h_amino_acids[AA_TRP].alt_atom_order[2] = ATOM_C;
    h_amino_acids[AA_TRP].alt_atom_order[3] = ATOM_CB;
    h_amino_acids[AA_TRP].alt_atom_order[4] = ATOM_O;
    h_amino_acids[AA_TRP].alt_atom_order[5] = ATOM_CG;
    h_amino_acids[AA_TRP].alt_atom_order[6] = ATOM_CD1;
    h_amino_acids[AA_TRP].alt_atom_order[7] = ATOM_CD2;
    h_amino_acids[AA_TRP].alt_atom_order[8] = ATOM_CE2;
    h_amino_acids[AA_TRP].alt_atom_order[9] = ATOM_CE3;
    h_amino_acids[AA_TRP].alt_atom_order[10] = ATOM_NE1;
    h_amino_acids[AA_TRP].alt_atom_order[11] = ATOM_CH2;
    h_amino_acids[AA_TRP].alt_atom_order[12] = ATOM_CZ2;
    h_amino_acids[AA_TRP].alt_atom_order[13] = ATOM_CZ3;
    SET_SC_ATOM(AA_TRP, 0, ATOM_O, ATOM_N, ATOM_CA, ATOM_C, 1.23f, 120.178f);
    SET_SC_ATOM(AA_TRP, 1, ATOM_CB, ATOM_O, ATOM_C, ATOM_CA, 1.53f, 110.852f);
    SET_SC_ATOM(AA_TRP, 2, ATOM_CG, ATOM_N, ATOM_CA, ATOM_CB, 1.50f, 114.10f);
    SET_SC_ATOM(AA_TRP, 3, ATOM_CD1, ATOM_CA, ATOM_CB, ATOM_CG, 1.36f, 126.712f);
    SET_SC_ATOM(AA_TRP, 4, ATOM_CD2, ATOM_CA, ATOM_CB, ATOM_CG, 1.44f, 126.712f);
    SET_SC_ATOM(AA_TRP, 5, ATOM_NE1, ATOM_CB, ATOM_CG, ATOM_CD1, 1.38f, 109.959f);
    SET_SC_ATOM(AA_TRP, 6, ATOM_CE2, ATOM_CB, ATOM_CG, ATOM_CD2, 1.41f, 107.842f);
    SET_SC_ATOM(AA_TRP, 7, ATOM_CE3, ATOM_CB, ATOM_CG, ATOM_CD2, 1.40f, 133.975f);
    SET_SC_ATOM(AA_TRP, 8, ATOM_CZ2, ATOM_CG, ATOM_CD2, ATOM_CE2, 1.40f, 120.0f);
    SET_SC_ATOM(AA_TRP, 9, ATOM_CZ3, ATOM_CG, ATOM_CD2, ATOM_CE3, 1.384f, 120.0f);
    SET_SC_ATOM(AA_TRP, 10, ATOM_CH2, ATOM_CD2, ATOM_CE2, ATOM_CZ2, 1.367f, 120.0f);

    // TYR - Tyrosine
    h_amino_acids[AA_TYR].n_atoms = 12;
    h_amino_acids[AA_TYR].n_sidechain_atoms = 9;
    h_amino_acids[AA_TYR].atom_order[0] = ATOM_N;
    h_amino_acids[AA_TYR].atom_order[1] = ATOM_CA;
    h_amino_acids[AA_TYR].atom_order[2] = ATOM_C;
    h_amino_acids[AA_TYR].atom_order[3] = ATOM_O;
    h_amino_acids[AA_TYR].atom_order[4] = ATOM_CB;
    h_amino_acids[AA_TYR].atom_order[5] = ATOM_CG;
    h_amino_acids[AA_TYR].atom_order[6] = ATOM_CD1;
    h_amino_acids[AA_TYR].atom_order[7] = ATOM_CD2;
    h_amino_acids[AA_TYR].atom_order[8] = ATOM_CE1;
    h_amino_acids[AA_TYR].atom_order[9] = ATOM_CE2;
    h_amino_acids[AA_TYR].atom_order[10] = ATOM_CZ;
    h_amino_acids[AA_TYR].atom_order[11] = ATOM_OH;
    h_amino_acids[AA_TYR].alt_atom_order[0] = ATOM_N;
    h_amino_acids[AA_TYR].alt_atom_order[1] = ATOM_CA;
    h_amino_acids[AA_TYR].alt_atom_order[2] = ATOM_C;
    h_amino_acids[AA_TYR].alt_atom_order[3] = ATOM_CB;
    h_amino_acids[AA_TYR].alt_atom_order[4] = ATOM_O;
    h_amino_acids[AA_TYR].alt_atom_order[5] = ATOM_CG;
    h_amino_acids[AA_TYR].alt_atom_order[6] = ATOM_CD1;
    h_amino_acids[AA_TYR].alt_atom_order[7] = ATOM_CD2;
    h_amino_acids[AA_TYR].alt_atom_order[8] = ATOM_CE1;
    h_amino_acids[AA_TYR].alt_atom_order[9] = ATOM_CE2;
    h_amino_acids[AA_TYR].alt_atom_order[10] = ATOM_OH;
    h_amino_acids[AA_TYR].alt_atom_order[11] = ATOM_CZ;
    SET_SC_ATOM(AA_TYR, 0, ATOM_O, ATOM_N, ATOM_CA, ATOM_C, 1.235f, 120.608f);
    SET_SC_ATOM(AA_TYR, 1, ATOM_CB, ATOM_O, ATOM_C, ATOM_CA, 1.53f, 110.852f);
    SET_SC_ATOM(AA_TYR, 2, ATOM_CG, ATOM_N, ATOM_CA, ATOM_CB, 1.51f, 113.744f);
    SET_SC_ATOM(AA_TYR, 3, ATOM_CD1, ATOM_CA, ATOM_CB, ATOM_CG, 1.39f, 120.937f);
    SET_SC_ATOM(AA_TYR, 4, ATOM_CD2, ATOM_CA, ATOM_CB, ATOM_CG, 1.39f, 120.937f);
    SET_SC_ATOM(AA_TYR, 5, ATOM_CE1, ATOM_CB, ATOM_CG, ATOM_CD1, 1.38f, 120.0f);
    SET_SC_ATOM(AA_TYR, 6, ATOM_CE2, ATOM_CB, ATOM_CG, ATOM_CD2, 1.38f, 120.0f);
    SET_SC_ATOM(AA_TYR, 7, ATOM_CZ, ATOM_CG, ATOM_CD1, ATOM_CE1, 1.378f, 120.0f);
    SET_SC_ATOM(AA_TYR, 8, ATOM_OH, ATOM_CD1, ATOM_CE1, ATOM_CZ, 1.375f, 120.0f);

    // VAL - Valine
    h_amino_acids[AA_VAL].n_atoms = 7;
    h_amino_acids[AA_VAL].n_sidechain_atoms = 4;
    h_amino_acids[AA_VAL].atom_order[0] = ATOM_N;
    h_amino_acids[AA_VAL].atom_order[1] = ATOM_CA;
    h_amino_acids[AA_VAL].atom_order[2] = ATOM_C;
    h_amino_acids[AA_VAL].atom_order[3] = ATOM_O;
    h_amino_acids[AA_VAL].atom_order[4] = ATOM_CB;
    h_amino_acids[AA_VAL].atom_order[5] = ATOM_CG1;
    h_amino_acids[AA_VAL].atom_order[6] = ATOM_CG2;
    h_amino_acids[AA_VAL].alt_atom_order[0] = ATOM_N;
    h_amino_acids[AA_VAL].alt_atom_order[1] = ATOM_CA;
    h_amino_acids[AA_VAL].alt_atom_order[2] = ATOM_C;
    h_amino_acids[AA_VAL].alt_atom_order[3] = ATOM_CB;
    h_amino_acids[AA_VAL].alt_atom_order[4] = ATOM_O;
    h_amino_acids[AA_VAL].alt_atom_order[5] = ATOM_CG1;
    h_amino_acids[AA_VAL].alt_atom_order[6] = ATOM_CG2;
    SET_SC_ATOM(AA_VAL, 0, ATOM_O, ATOM_N, ATOM_CA, ATOM_C, 1.235f, 120.472f);
    SET_SC_ATOM(AA_VAL, 1, ATOM_CB, ATOM_O, ATOM_C, ATOM_CA, 1.54f, 111.381f);
    SET_SC_ATOM(AA_VAL, 2, ATOM_CG1, ATOM_N, ATOM_CA, ATOM_CB, 1.52f, 110.7f);
    SET_SC_ATOM(AA_VAL, 3, ATOM_CG2, ATOM_N, ATOM_CA, ATOM_CB, 1.52f, 110.4f);

    // UNK - Unknown. The CPU reference (AminoAcid('X', "UNK", "Unknown") in
    // amino_acid.h) has no sidechain map entries, so it emits only the 3
    // backbone atoms (no O). Do NOT alias to GLY here: GLY's O is placed via
    // the same sidechain-geometry mechanism as a real sidechain atom, so
    // aliasing would give UNK a 4th atom the CPU never writes.
    h_amino_acids[AA_UNK].n_atoms = 3;
    h_amino_acids[AA_UNK].n_sidechain_atoms = 0;
    h_amino_acids[AA_UNK].atom_order[0] = ATOM_N;
    h_amino_acids[AA_UNK].atom_order[1] = ATOM_CA;
    h_amino_acids[AA_UNK].atom_order[2] = ATOM_C;
    h_amino_acids[AA_UNK].alt_atom_order[0] = ATOM_N;
    h_amino_acids[AA_UNK].alt_atom_order[1] = ATOM_CA;
    h_amino_acids[AA_UNK].alt_atom_order[2] = ATOM_C;

#undef SET_SC_ATOM
}

void initGPUAminoAcidTable() {
    int device = 0;
    cudaGetDevice(&device);
    std::lock_guard<std::mutex> lock(g_aa_table_mutex);
    if (g_aa_table_initialized_devices.count(device)) return;

    GPUAminoAcid h_amino_acids[GPU_NUM_AMINO_ACIDS + 1];
    buildHostAminoAcidTable(h_amino_acids);

    cudaError_t err = cudaMemcpyToSymbol(d_amino_acids, h_amino_acids,
        sizeof(GPUAminoAcid) * (GPU_NUM_AMINO_ACIDS + 1));
    if (err != cudaSuccess) {
        throw std::runtime_error(std::string("Failed to copy amino acid table to GPU: ") +
                                 cudaGetErrorString(err));
    }

    // Initialize atom name string table
    char h_atom_name_strs[36][5] = {};
    for (int i = 0; i < 36; i++) {
        strncpy(h_atom_name_strs[i], ATOM_NAME_STRINGS[i], 4);
        h_atom_name_strs[i][4] = '\0';
    }
    cudaMemcpyToSymbol(d_atom_name_strs, h_atom_name_strs, sizeof(h_atom_name_strs));

    // Initialize residue name string table (indexed by AminoAcidIndex)
    static const char* RESIDUE_NAME_STRS[21] = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "UNK"};
    char h_residue_name_strs[21][4] = {};
    for (int i = 0; i < 21; i++) {
        strncpy(h_residue_name_strs[i], RESIDUE_NAME_STRS[i], 3);
        h_residue_name_strs[i][3] = '\0';
    }
    cudaMemcpyToSymbol(d_residue_name_strs, h_residue_name_strs, sizeof(h_residue_name_strs));

    // Initialize per-AA scalar lookup tables for buildSidechainMeta kernels
    // Bucket ID formula: sc<=2→0, sc<=4→1, sc<=6→2, else→3 (matches AA_BUCKET_ID exactly)
    {
        int8_t h_sct[21], h_ac[21], h_bid[21];
        for (int i = 0; i < GPU_NUM_AMINO_ACIDS; i++) {
            int sc = SIDECHAIN_TORSION_COUNTS[i];
            h_sct[i] = (int8_t)sc;
            h_ac[i] = (int8_t)ATOM_COUNTS[i];
            h_bid[i] = (sc <= 2) ? 0 : (sc <= 4) ? 1
                                   : (sc <= 6)   ? 2
                                                 : 3;
        }
        // AA_UNK (index 20): 3 backbone-only atoms, no sidechain torsions —
        // matches the h_amino_acids[AA_UNK] entry above, not GLY (see there).
        h_sct[AA_UNK] = 0;
        h_ac[AA_UNK] = 3;
        h_bid[AA_UNK] = 0;
        cudaMemcpyToSymbol(d_sc_torsion_counts_k, h_sct, 21 * sizeof(int8_t));
        cudaMemcpyToSymbol(d_atom_counts_k, h_ac, 21 * sizeof(int8_t));
        cudaMemcpyToSymbol(d_bucket_id_k, h_bid, 21 * sizeof(int8_t));
    }

    g_aa_table_initialized_devices.insert(device);
}

void cleanupGPUAminoAcidTable() {
    std::lock_guard<std::mutex> lock(g_aa_table_mutex);
    g_aa_table_initialized_devices.clear();
}

bool isGPUAminoAcidTableInitialized() {
    int device = 0;
    cudaGetDevice(&device);
    std::lock_guard<std::mutex> lock(g_aa_table_mutex);
    return g_aa_table_initialized_devices.count(device) != 0;
}

// ============================================================================
// Device Helper Functions
// ============================================================================

__device__ inline float3 cross_product_dev(const float3& a, const float3& b) {
    return make_float3(
        a.y * b.z - b.y * a.z,
        a.z * b.x - b.z * a.x,
        a.x * b.y - b.x * a.y);
}

__device__ inline float norm_dev(const float3& v) {
    return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}

/**
 * @brief Place a single atom using NeRF algorithm (device function)
 */
__device__ float3 place_atom_dev(
    const float3& a, // Atom A
    const float3& b, // Atom B
    const float3& c, // Atom C
    float bond_length,
    float bond_angle,   // in degrees
    float torsion_angle // in degrees
) {
    // Convert angles to radians
    float ba_rad = bond_angle * M_PI / 180.0f;
    float ta_rad = torsion_angle * M_PI / 180.0f;

    // Calculate vectors AB and BC
    float3 ab = make_float3(b.x - a.x, b.y - a.y, b.z - a.z);
    float3 bc = make_float3(c.x - b.x, c.y - b.y, c.z - b.z);

    // Normalize BC
    float bc_norm = norm_dev(bc);
    float3 bcn = make_float3(bc.x / bc_norm, bc.y / bc_norm, bc.z / bc_norm);

    // Calculate current atom position in local frame
    float curr_x = -bond_length * cosf(ba_rad);
    float curr_y = bond_length * cosf(ta_rad) * sinf(ba_rad);
    float curr_z = bond_length * sinf(ta_rad) * sinf(ba_rad);

    // Calculate cross product: n = AB x BCn
    float3 n = cross_product_dev(ab, bcn);

    // Normalize n
    float n_norm = norm_dev(n);
    n.x /= n_norm;
    n.y /= n_norm;
    n.z /= n_norm;

    // Calculate nbc = n x BCn
    float3 nbc = cross_product_dev(n, bcn);

    // Build rotation matrix M = [BCn, nbc, n]
    // Apply rotation and translation
    return make_float3(
        bcn.x * curr_x + nbc.x * curr_y + n.x * curr_z + c.x,
        bcn.y * curr_x + nbc.y * curr_y + n.y * curr_z + c.y,
        bcn.z * curr_x + nbc.z * curr_y + n.z * curr_z + c.z);
}

/**
 * @brief Get atom coordinates by name index from reconstructed atoms
 */
__device__ void get_atom_coords_by_name(
    const float* atom_coords,    // All atom coordinates for this residue
    const int8_t* placed_atoms,  // Which atoms have been placed
    int n_placed,                // Number of atoms placed so far
    int8_t target_atom,          // Atom name index to find
    const GPUAminoAcid& aa,      // Amino acid info
    float& x, float& y, float& z // Output coordinates
) {
    // First check backbone atoms (always in positions 0-2)
    if (target_atom == ATOM_N) {
        x = atom_coords[0];
        y = atom_coords[1];
        z = atom_coords[2];
        return;
    }
    if (target_atom == ATOM_CA) {
        x = atom_coords[3];
        y = atom_coords[4];
        z = atom_coords[5];
        return;
    }
    if (target_atom == ATOM_C) {
        x = atom_coords[6];
        y = atom_coords[7];
        z = atom_coords[8];
        return;
    }

    // Search in placed sidechain atoms
    for (int i = 0; i < n_placed; i++) {
        if (placed_atoms[i] == target_atom) {
            // Sidechain atoms start at index 3 (after N, CA, C)
            int idx = (3 + i) * 3;
            x = atom_coords[idx];
            y = atom_coords[idx + 1];
            z = atom_coords[idx + 2];
            return;
        }
    }

    // Fallback: should not reach here
    x = y = z = 0.0f;
}

// ============================================================================
// Bucket classification for reduced warp divergence
// ============================================================================

// Bucket classification: sc<=2→0, sc<=4→1, sc<=6→2, else→3
// (initialized into d_bucket_id_k constant memory in initGPUAminoAcidTable)

// ============================================================================
// Device helper: inline backbone + alt-order lookup (replaces d_bb_full scatter)
// ============================================================================

/**
 * Given global residue index res_idx, find its struct si via binary search over d_rfo,
 * then return the backbone pointer (prevAtoms for r_local==0, backbone output otherwise)
 * and the per-struct alt_order flag.
 */
__device__ inline void get_bb_and_alt(
    int res_idx, int n_structs,
    const float* __restrict__ d_bb_out,
    const float* __restrict__ d_pa,
    const int* __restrict__ d_rfo,
    const int* __restrict__ d_bbo,
    const int8_t* __restrict__ d_use_alt_per_struct,
    const float*& bb_ptr,
    int& alt_order) {
    int lo = 0, hi = n_structs - 1;
    while (lo < hi) {
        int mid = (lo + hi + 1) >> 1;
        if (d_rfo[mid] <= res_idx)
            lo = mid;
        else
            hi = mid - 1;
    }
    int si = lo;
    int r_local = res_idx - d_rfo[si];
    bb_ptr = (r_local == 0) ? (d_pa + si * 9)
                            : (d_bb_out + d_bbo[si] + (r_local - 1) * 9);
    alt_order = (int)d_use_alt_per_struct[si];
}

// ============================================================================
// CUDA Kernel for Sidechain Reconstruction
// ============================================================================

// ============================================================================
// Templated bucket kernels: MAX_SC reduces warp divergence within each bucket.
// ============================================================================

/**
 * @brief Sidechain reconstruction kernel for a single bucket of residues.
 *
 * MAX_SC is the maximum n_sidechain_atoms for any AA in this bucket.
 * Compile-time sizing of atom_coords[] and placed_atoms[] allows the compiler
 * to keep them in registers and unroll the sidechain loop.
 *
 * @param bucket_id  Bucket selector in the range [0, 3]
 * @param d_bsz_all  Per-bucket residue counts
 * @param d_boff_all Per-bucket offsets into d_bidx
 * @param d_bidx     Global residue indices grouped by bucket
 */
template <int MAX_SC>
__global__ void reconstruct_sidechains_kernel(
    int bucket_id,
    const int* __restrict__ d_bsz_all,
    const int* __restrict__ d_boff_all,
    const int* __restrict__ d_bidx,
    int n_structs,
    const int8_t* __restrict__ residue_types,
    const float* __restrict__ d_bb_out, // backbone output (residues 1..n-1)
    const float* __restrict__ d_pa,     // prevAtoms per struct [n_structs*9]
    const int* __restrict__ d_rfo,      // residue full offsets [n_structs+1]
    const int* __restrict__ d_bbo,      // bb float offsets [n_structs+1]
    const float* __restrict__ sidechain_angles,
    const int2* __restrict__ d_offsets,            // [n_residues+1] .x=sc angle off, .y=atom off
    const int8_t* __restrict__ use_alt_per_struct, // [n_structs]
    float* __restrict__ output_coords,
    int8_t* __restrict__ output_atom_names) {
    int n_bucket = d_bsz_all[bucket_id];
    int bid = blockIdx.x * blockDim.x + threadIdx.x;
    if (bid >= n_bucket) return;
    int res_idx = d_bidx[d_boff_all[bucket_id] + bid];

    int8_t aa_type = residue_types[res_idx];
    if (aa_type < 0 || aa_type > GPU_NUM_AMINO_ACIDS) aa_type = AA_UNK;
    const GPUAminoAcid& aa = d_amino_acids[aa_type];

    int out_offset = d_offsets[res_idx].y;
    float* my_output = output_coords + out_offset * 3;
    int8_t* my_atom_names = output_atom_names + out_offset;

    const float* bb;
    int alt_order;
    get_bb_and_alt(res_idx, n_structs, d_bb_out, d_pa, d_rfo, d_bbo,
        use_alt_per_struct, bb, alt_order);

    float atom_coords[(3 + MAX_SC) * 3];
    atom_coords[0] = bb[0];
    atom_coords[1] = bb[1];
    atom_coords[2] = bb[2];
    atom_coords[3] = bb[3];
    atom_coords[4] = bb[4];
    atom_coords[5] = bb[5];
    atom_coords[6] = bb[6];
    atom_coords[7] = bb[7];
    atom_coords[8] = bb[8];

    int8_t placed_atoms[MAX_SC];
    int n_placed = 0;
    int sc_start = d_offsets[res_idx].x;

    for (int i = 0; i < aa.n_sidechain_atoms; i++) {
        const GPUSidechainAtom& sc_atom = aa.sidechain[i];
        float ax, ay, az, bx, by, bz, cx, cy, cz;
        get_atom_coords_by_name(atom_coords, placed_atoms, n_placed,
            sc_atom.dep_atoms[0], aa, ax, ay, az);
        get_atom_coords_by_name(atom_coords, placed_atoms, n_placed,
            sc_atom.dep_atoms[1], aa, bx, by, bz);
        get_atom_coords_by_name(atom_coords, placed_atoms, n_placed,
            sc_atom.dep_atoms[2], aa, cx, cy, cz);
        float torsion = sidechain_angles[sc_start + i];
        float3 d = place_atom_dev(
            make_float3(ax, ay, az), make_float3(bx, by, bz), make_float3(cx, cy, cz),
            sc_atom.bond_length, sc_atom.bond_angle, torsion);
        int idx = (3 + n_placed) * 3;
        atom_coords[idx] = d.x;
        atom_coords[idx + 1] = d.y;
        atom_coords[idx + 2] = d.z;
        placed_atoms[n_placed++] = sc_atom.atom_name;
    }

    const int8_t* atom_order = (alt_order == 1) ? aa.alt_atom_order : aa.atom_order;

    for (int i = 0; i < aa.n_atoms; i++) {
        int8_t atom_name = atom_order[i];
        float x, y, z;
        if (atom_name == ATOM_N) {
            x = atom_coords[0];
            y = atom_coords[1];
            z = atom_coords[2];
        } else if (atom_name == ATOM_CA) {
            x = atom_coords[3];
            y = atom_coords[4];
            z = atom_coords[5];
        } else if (atom_name == ATOM_C) {
            x = atom_coords[6];
            y = atom_coords[7];
            z = atom_coords[8];
        } else {
            bool found = false;
            for (int j = 0; j < n_placed && !found; j++) {
                if (placed_atoms[j] == atom_name) {
                    int idx = (3 + j) * 3;
                    x = atom_coords[idx];
                    y = atom_coords[idx + 1];
                    z = atom_coords[idx + 2];
                    found = true;
                }
            }
            if (!found) x = y = z = 0.0f;
        }
        my_output[i * 3] = x;
        my_output[i * 3 + 1] = y;
        my_output[i * 3 + 2] = z;
        my_atom_names[i] = atom_name;
    }
}

// Explicit instantiations for the four bucket sizes.
template __global__ void reconstruct_sidechains_kernel<2>(
    int, const int*, const int*, const int*, int, const int8_t*, const float*, const float*,
    const int*, const int*, const float*, const int2*, const int8_t*, float*, int8_t*);
template __global__ void reconstruct_sidechains_kernel<4>(
    int, const int*, const int*, const int*, int, const int8_t*, const float*, const float*,
    const int*, const int*, const float*, const int2*, const int8_t*, float*, int8_t*);
template __global__ void reconstruct_sidechains_kernel<6>(
    int, const int*, const int*, const int*, int, const int8_t*, const float*, const float*,
    const int*, const int*, const float*, const int2*, const int8_t*, float*, int8_t*);
template __global__ void reconstruct_sidechains_kernel<11>(
    int, const int*, const int*, const int*, int, const int8_t*, const float*, const float*,
    const int*, const int*, const float*, const int2*, const int8_t*, float*, int8_t*);

// ============================================================================
// AtomCoordinate Fill Kernel
// ============================================================================

/**
 * One thread per residue. Each thread writes all atoms for that residue into the
 * output AtomCoordinate array. Atom names and residue names come from __constant__
 * lookup tables, avoiding any host-pointer access in device code.
 *
 * output[offsets[res]] .. output[offsets[res+1]-1] are the atoms for residue res.
 */
__global__ void fill_atom_coordinates_kernel(
    int n_residues,
    int n_structs,
    int n_tf,
    const int8_t* __restrict__ residue_types,      // [n_residues] AminoAcidIndex
    const float* __restrict__ d_tf,                // [n_tf] device TF floats (indexed by global res_idx)
    const int* __restrict__ struct_res_offsets,    // [n_structs + 1] cumulative residue counts
    const int* __restrict__ base_anum_per_struct,  // [n_structs] header.idxAtom per structure
    const int* __restrict__ base_ridx_per_struct,  // [n_structs] header.idxResidue per structure
    const int* __restrict__ base_delta_per_struct, // [n_structs] ac_offset - gpu_off_start per structure
    const char* __restrict__ chain_per_struct,     // [n_structs * CHAIN_ID_LENGTH] chain bytes per structure (mirrors ChainId/FixedStr<CHAIN_ID_LENGTH>)
    const int2* __restrict__ offsets,              // [n_residues + 1] .y = atom prefix-sum
    const int8_t* __restrict__ atom_names,         // [total_atoms] AtomNameIndex
    const float* __restrict__ coords,              // [total_atoms * 3]
    AtomCoordinate* __restrict__ output            // [total_atoms + n_oxt], caller-allocated
) {
    int res_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (res_idx >= n_residues) return;

    int8_t aa = residue_types[res_idx];
    if (aa < 0 || aa > AA_UNK) aa = AA_GLY;

    // Binary search to find which structure this residue belongs to
    int lo = 0, hi = n_structs - 1;
    while (lo < hi) {
        int mid = (lo + hi) >> 1;
        if (struct_res_offsets[mid + 1] <= res_idx)
            lo = mid + 1;
        else
            hi = mid;
    }
    int si = lo;

    int local_idx = res_idx - struct_res_offsets[si];
    float tf = (res_idx < n_tf && d_tf) ? d_tf[res_idx] : 0.0f;
    int ri = base_ridx_per_struct[si] + local_idx;
    int gpu_off_start = offsets[struct_res_offsets[si]].y;
    int anum = base_anum_per_struct[si] + offsets[res_idx].y - gpu_off_start;
    const char* ch4 = chain_per_struct + si * CHAIN_ID_LENGTH;
    int out_base = base_delta_per_struct[si]; // shifts into correct structure slot

    int out_start = offsets[res_idx].y;
    int out_end = offsets[res_idx + 1].y;

    for (int a = out_start; a < out_end; a++) {
        int8_t name_idx = atom_names[a];
        if (name_idx < 0 || name_idx >= 36) name_idx = 0;

        AtomCoordinate& out = output[out_base + a];

        // atom name (FixedStr<4> = char[5])
        const char* an = d_atom_name_strs[name_idx];
        out.atom.data[0] = an[0];
        out.atom.data[1] = an[1];
        out.atom.data[2] = an[2];
        out.atom.data[3] = an[3];
        out.atom.data[4] = '\0';

        // residue name (FixedStr<3> = char[4])
        const char* rn = d_residue_name_strs[aa];
        out.residue.data[0] = rn[0];
        out.residue.data[1] = rn[1];
        out.residue.data[2] = rn[2];
        out.residue.data[3] = '\0';

        // chain (FixedStr<CHAIN_ID_LENGTH> = char[CHAIN_ID_LENGTH+1])
        for (int cb = 0; cb < CHAIN_ID_LENGTH; cb++) out.chain.data[cb] = ch4[cb];
        out.chain.data[CHAIN_ID_LENGTH] = '\0';

        out.atom_index = anum + (a - out_start);
        out.residue_index = ri;
        out.coordinate.x = coords[a * 3 + 0];
        out.coordinate.y = coords[a * 3 + 1];
        out.coordinate.z = coords[a * 3 + 2];
        out.occupancy = 0.0f;
        out.tempFactor = tf;
        out.model = 1;
        out.insertion_code = ' ';
        out.altloc = ' ';
    }
}

// ============================================================================
// buildSidechainMeta kernels
// ============================================================================

// Kernel: compute per-residue SC torsion count and atom count from d_rt.
// We launch n_residues+1 threads; thread n writes sentinel zeros so that
// CUB ExclusiveScan(n+1) produces the full n+1 offset array in one call.
// d_meta[i].x = sc torsion count, d_meta[i].y = atom count.
__global__ void sc_meta_counts_kernel(
    const int8_t* d_rt, int n_residues,
    int2* d_meta) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i > n_residues) return;
    if (i == n_residues) {
        d_meta[i] = {0, 0};
        return;
    }
    int8_t aa = d_rt[i];
    // Valid range is [0, GPU_NUM_AMINO_ACIDS] inclusive: AA_UNK == GPU_NUM_AMINO_ACIDS
    // (20) is itself a valid index with its own d_atom_counts_k/d_sc_torsion_counts_k
    // row (3 atoms, 0 torsions — see buildHostAminoAcidTable). A `>=` here would
    // misroute AA_UNK to AA_GLY's row (4 atoms), corrupting every UNK residue's
    // reconstructed atom count.
    if (aa < 0 || aa > (int8_t)GPU_NUM_AMINO_ACIDS) aa = AA_GLY;
    d_meta[i] = {(int)d_sc_torsion_counts_k[aa], (int)d_atom_counts_k[aa]};
}

// Kernel: count residues per bucket using atomicAdd into d_bsz[4].
// d_bsz must be zeroed before launch.
__global__ void bucket_count_kernel(
    const int8_t* d_rt, int n_residues, int* d_bsz) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_residues) return;
    int8_t aa = d_rt[i];
    if (aa < 0 || aa > (int8_t)GPU_NUM_AMINO_ACIDS) aa = AA_GLY;
    atomicAdd(&d_bsz[(int)d_bucket_id_k[aa]], 1);
}

// Kernel: scatter residue indices into per-bucket arrays.
// d_boff[4] = base offset per bucket into d_bucket_out.
// d_bcnt[4] must be zeroed before launch; used as per-bucket atomic position counter.
__global__ void bucket_scatter_kernel(
    const int8_t* d_rt, int n_residues,
    int* d_bucket_out, const int* d_boff, int* d_bcnt) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_residues) return;
    int8_t aa = d_rt[i];
    if (aa < 0 || aa > (int8_t)GPU_NUM_AMINO_ACIDS) aa = AA_GLY;
    int bid = (int)d_bucket_id_k[aa];
    int pos = d_boff[bid] + atomicAdd(&d_bcnt[bid], 1);
    d_bucket_out[pos] = i;
}

// Kernel: compute d_boff[4] as exclusive prefix sum of d_bsz[4], then zero d_bcnt[4].
// Runs as a single thread — 4 additions, negligible cost, eliminates D→H round-trip.
__global__ void bucket_prefix_sum_kernel(const int* d_bsz, int* d_boff, int* d_bcnt) {
    d_boff[0] = 0;
    d_boff[1] = d_bsz[0];
    d_boff[2] = d_bsz[0] + d_bsz[1];
    d_boff[3] = d_bsz[0] + d_bsz[1] + d_bsz[2];
    d_bcnt[0] = d_bcnt[1] = d_bcnt[2] = d_bcnt[3] = 0;
}

// ============================================================================
// buildSidechainMeta
// ============================================================================

// Binary op for CUB ExclusiveScan over int2 (adds .x and .y independently).
struct Int2Add {
    __device__ __host__ int2 operator()(int2 a, int2 b) const {
        return {a.x + b.x, a.y + b.y};
    }
};

SCMetaPtrs buildSidechainMeta(SidechainBuffers& ctx, const int8_t* d_rt, int n_residues,
    cudaStream_t stream) {
    if (!isGPUAminoAcidTableInitialized()) initGPUAminoAcidTable();
    const int n1 = n_residues + 1; // +1 for sentinel → ExclusiveScan gives n+1 outputs

    int2* d_meta_counts = ctx.d_meta_counts.ensure_n<int2>(n1, stream);
    int2* d_meta_off = ctx.d_meta_off.ensure_n<int2>(n1, stream);

    // Step 1: compute per-residue counts — .x=sc torsion count, .y=atom count
    const int BS = 256;
    int blocks = (n1 + BS - 1) / BS;
    sc_meta_counts_kernel<<<blocks, BS, 0, stream>>>(d_rt, n_residues, d_meta_counts);

    // Step 2: single CUB ExclusiveScan over int2 — one kernel launch instead of two.
    // First call (d_temp=nullptr) queries required temp storage size (no GPU work).
    {
        size_t temp_sz = 0;
        cub::DeviceScan::ExclusiveScan(nullptr, temp_sz,
            d_meta_counts, d_meta_off, Int2Add{}, int2{0, 0}, n1, stream);
        void* d_temp = ctx.d_cub_temp.ensure(temp_sz ? temp_sz : 1, stream);
        cub::DeviceScan::ExclusiveScan(d_temp, temp_sz,
            d_meta_counts, d_meta_off, Int2Add{}, int2{0, 0}, n1, stream);
    }

    SCMetaPtrs meta;
    meta.d_offsets = d_meta_off;
    // Layout in d_bsz_boff_bcnt: [bsz0..3 | boff0..3 | bcnt0..3]  (12 ints)
    int* d_bsz = ctx.d_bsz_boff_bcnt.ensure_n<int>(12, stream);
    int* d_boff = d_bsz + 4;
    int* d_bcnt = d_boff + 4;
    int* d_bidx = ctx.d_bucket_indices.ensure_n<int>(n_residues, stream);

    // Zero bsz before counting (bcnt is zeroed by the prefix-sum kernel).
    cudaMemsetAsync(d_bsz, 0, 4 * sizeof(int), stream);

    // Step 3: count residues per bucket.
    bucket_count_kernel<<<(n_residues + BS - 1) / BS, BS, 0, stream>>>(d_rt, n_residues, d_bsz);

    // Step 4: compute d_boff[4] as the exclusive prefix sum of d_bsz[4] entirely on device,
    //         and zero d_bcnt[4].
    bucket_prefix_sum_kernel<<<1, 1, 0, stream>>>(d_bsz, d_boff, d_bcnt);

    // Step 5: scatter residue indices into bucket arrays.
    bucket_scatter_kernel<<<(n_residues + BS - 1) / BS, BS, 0, stream>>>(
        d_rt, n_residues, d_bidx, d_boff, d_bcnt);

    meta.d_bsz = d_bsz;
    meta.d_boff = d_boff;
    meta.d_bidx = d_bidx;

    return meta;
}

// ============================================================================
// Host Wrapper Functions
// ============================================================================

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
    cudaStream_t stream) {
    if (n_residues == 0) return;

    // meta.d_output_offsets[n_residues] holds total atom count; it was computed by buildSidechainMeta
    // on device. We don't download it — output_coords is sized from the caller's CPU value.
    // Compute a safe upper bound: n_residues * GPU_MAX_ATOMS_PER_AA.
    const int total_atoms_ub = n_residues * GPU_MAX_ATOMS_PER_AA;

    auto* d_output_coords = ctx.output_coords.ensure_n<float>(total_atoms_ub * 3, stream);
    auto* d_output_atom_names = ctx.output_atom_names.ensure_n<int8_t>(total_atoms_ub, stream);
    constexpr int BLOCK_SIZE = 256;

    // Upload use_alt_per_struct (≤100 bytes) — the only remaining H→D transfer here
    auto* h_alt = ctx.ph_alt.ensure_n<int8_t>(n_structs);
    auto* d_alt = ctx.d_alt.ensure_n<int8_t>(n_structs, stream);
    memcpy(h_alt, use_alt_per_struct.data(), (size_t)n_structs);
    cudaMemcpyAsync(d_alt, h_alt, (size_t)n_structs, cudaMemcpyHostToDevice, stream);

    // Cache device pointers for fillAtomCoordinatesGPU (called next on same stream)
    ctx.d_rt_cached = d_rt;
    ctx.d_ooff_cached = meta.d_offsets;

    // Upper-bound grid: each kernel reads its n_bucket from device and exits early if tid >= n_bucket.
    const int nb_ub = (n_residues + BLOCK_SIZE - 1) / BLOCK_SIZE;
#define LAUNCH_BUCKET(MAX_SC, BID) \
    reconstruct_sidechains_kernel<MAX_SC><<<nb_ub, BLOCK_SIZE, 0, stream>>>( \
        BID, meta.d_bsz, meta.d_boff, meta.d_bidx, n_structs, \
        d_rt, d_bb_out, d_pa, d_rfo, d_bbo, \
        d_sidechain_angles, meta.d_offsets, d_alt, \
        d_output_coords, d_output_atom_names);

        LAUNCH_BUCKET(2, 0)  // GLY, ALA
        LAUNCH_BUCKET(4, 1)  // SER, CYS, PRO, THR, VAL
        LAUNCH_BUCKET(6, 2)  // ASN, ASP, ILE, LEU, MET, GLN, GLU, LYS
        LAUNCH_BUCKET(11, 3) // HIS, PHE, ARG, TYR, TRP

#undef LAUNCH_BUCKET
    // No sync: fill kernel on the same stream will wait automatically.
}

// ============================================================================
// fillAtomCoordinatesGPU
// ============================================================================
/**
 * Fills an AtomCoordinate array on the GPU using the outputs that remain on-device
 * from the preceding reconstructSidechainsGPU call (ctx.output_coords, etc.).
 * Must be called immediately after reconstructSidechainsGPU.
 *
 * Per-residue metadata (temp_factors, res_local_idx, res_atom_num_start, res_chain)
 * is uploaded from host. The filled AtomCoordinate array is downloaded into output[].
 *
 * total_atoms must equal output_offsets.back() (the value returned by reconstructSidechainsGPU).
 */
void fillAtomCoordinatesGPU(
    SidechainBuffers& ctx,
    int n_residues,
    int n_structs,
    int total_atoms,
    int total_output_atoms,                        // >= total_atoms (includes OXT slots)
    int n_tf,                                      // length of d_tf array
    const float* d_tf,                             // [n_tf] device TF ptr (from continuize_all)
    const std::vector<int>& struct_res_offsets,    // [n_structs + 1] cumulative residue counts
    const std::vector<int>& base_anum_per_struct,  // [n_structs] header.idxAtom per struct
    const std::vector<int>& base_ridx_per_struct,  // [n_structs] header.idxResidue per struct
    const std::vector<int>& base_delta_per_struct, // [n_structs] ac_offset - gpu_off_start per struct
    const std::vector<char>& chain_per_struct,     // [n_structs * CHAIN_ID_LENGTH] chain bytes per struct
    AtomCoordinate* output,                        // pre-allocated pinned [total_output_atoms]
    cudaStream_t stream_a,                         // compute stream for this pipeline
    cudaStream_t stream_b,                         // D→H stream for this pipeline
    cudaEvent_t e_fill_done,                       // signals stream_b after fill kernel
    cudaEvent_t e_slot_ready,                      // guards next reuse of output buffer
    cudaEvent_t cpu_sync_event,                    // per-batch: CPU sync after D→H
    bool skip_download                             // run fill kernel but skip the D→H
) {
    if (n_residues == 0 || total_atoms == 0) return;

    // Pack [int[n+1] sro | int[n] banum | int[n] bridx | int[n] bdelta |
    //       char[n*4] chain] into one transfer.
    // Device output buffer is separate (large, kernel write target).

    const size_t sz_sro = (size_t)(n_structs + 1) * sizeof(int);
    const size_t sz_banum = (size_t)n_structs * sizeof(int);
    const size_t sz_bridx = (size_t)n_structs * sizeof(int);
    const size_t sz_bdelta = (size_t)n_structs * sizeof(int);
    const size_t sz_chain = (size_t)n_structs * CHAIN_ID_LENGTH * sizeof(char);
    const size_t off_sro = 0;
    const size_t off_banum = off_sro + sz_sro;
    const size_t off_bridx = off_banum + sz_banum;
    const size_t off_bdel = off_bridx + sz_bridx;
    const size_t off_chain = off_bdel + sz_bdelta; // char: 1-byte, no alignment gap
    const size_t fill_sz = off_chain + sz_chain;

    char* hf = static_cast<char*>(ctx.ph_fill_pack.ensure(fill_sz));
    char* df = static_cast<char*>(ctx.d_fill_pack.ensure(fill_sz, stream_a));

    int* h_sro = reinterpret_cast<int*>(hf + off_sro);
    int* h_banum = reinterpret_cast<int*>(hf + off_banum);
    int* h_bridx = reinterpret_cast<int*>(hf + off_bridx);
    int* h_bdelta = reinterpret_cast<int*>(hf + off_bdel);
    char* h_chain = reinterpret_cast<char*>(hf + off_chain);
    int* d_struct_res_off = reinterpret_cast<int*>(df + off_sro);
    int* d_base_anum = reinterpret_cast<int*>(df + off_banum);
    int* d_base_ridx = reinterpret_cast<int*>(df + off_bridx);
    int* d_base_delta = reinterpret_cast<int*>(df + off_bdel);
    char* d_chain_ps = reinterpret_cast<char*>(df + off_chain);

    // Guard: wait for any prior D→H using this output buffer to complete before overwriting it.
    cudaStreamWaitEvent(stream_a, e_slot_ready);

    auto* d_output = ctx.d_ac_output.ensure_n<AtomCoordinate>(total_output_atoms, stream_a);

    memcpy(h_sro, struct_res_offsets.data(), sz_sro);
    memcpy(h_banum, base_anum_per_struct.data(), sz_banum);
    memcpy(h_bridx, base_ridx_per_struct.data(), sz_bridx);
    memcpy(h_bdelta, base_delta_per_struct.data(), sz_bdelta);
    memcpy(h_chain, chain_per_struct.data(), sz_chain);

    // Single transfer for all five fill-kernel metadata arrays (~2 KB for 100 structs)
    cudaMemcpyAsync(df, hf, fill_sz, cudaMemcpyHostToDevice, stream_a);

    constexpr int BLOCK_SIZE = 256;
    int num_blocks = (n_residues + BLOCK_SIZE - 1) / BLOCK_SIZE;

    fill_atom_coordinates_kernel<<<num_blocks, BLOCK_SIZE, 0, stream_a>>>(
        n_residues,
        n_structs,
        n_tf,
        ctx.d_rt_cached, // set by reconstructSidechainsGPU
        d_tf,
        d_struct_res_off,
        d_base_anum,
        d_base_ridx,
        d_base_delta,
        d_chain_ps,
        ctx.d_ooff_cached, // set by reconstructSidechainsGPU
        ctx.output_atom_names.ensure_n<int8_t>(total_atoms, stream_a),
        ctx.output_coords.ensure_n<float>(total_atoms * 3, stream_a),
        d_output);

    // Signal stream A: fill kernel is done for this slot.
    cudaEventRecord(e_fill_done, stream_a);

    if (skip_download) {
        // Caller consumes ctx.d_ac_output on the device (e.g. GPU PDB writer).
        // Record slot_ready on stream A so the next reuse guard is satisfied.
        cudaEventRecord(e_slot_ready, stream_a);
        return;
    }

    // Stream B: wait for fill kernel, then DMA result to pinned host buffer.
    cudaStreamWaitEvent(stream_b, e_fill_done);
    cudaMemcpyAsync(output, d_output, total_output_atoms * sizeof(AtomCoordinate),
        cudaMemcpyDeviceToHost, stream_b);
    // Per-slot guard: unblocks stream A's next use of this output buffer slot.
    cudaEventRecord(e_slot_ready, stream_b);
    // Per-batch CPU-sync: caller's cudaEventSynchronize(cpu_sync_event) waits on this.
    cudaEventRecord(cpu_sync_event, stream_b);
    // Caller must cudaEventSynchronize(cpu_sync_event) before reading `output`.
}
