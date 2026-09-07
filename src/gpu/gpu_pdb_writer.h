// SPDX-License-Identifier: MIT
/**
 * File: gpu_pdb_writer.h
 * Project: foldcomp
 * Description:
 *     GPU PDB writer: format the per-atom PDB ATOM records on the device (the
 *     O(N_atoms) work that otherwise bottlenecks the host text path), leaving the
 *     rare structural records (title/MODEL/TER) to the host stitcher
 *     (writeSegmentsToPDB with its precomputed-lines parameter).
 */
#pragma once

#ifdef FOLDCOMP_WITH_CUDA

#include "structure_codec.h" // AtomCoordinate

// Forward-declare the CUDA stream handle instead of including <cuda_runtime.h>.
// cudaStream_t is an opaque pointer (typedef struct CUstream_st* cudaStream_t),
// so naming the exact type keeps this header CUDA-include-free while remaining
// type-safe. C++ permits this identical typedef to coexist with the one from
// <cuda_runtime.h> when a translation unit includes both.
typedef struct CUstream_st* cudaStream_t;

// In-pipeline path: format ATOM lines directly from a device AtomCoordinate array
// (no host round-trip). Writes n_atoms * PDB_ATOM_LINE_LEN bytes into `d_lines`
// (device) on `stream`. Caller owns/sizes both buffers.
void formatAtomLinesDeviceGPU(const AtomCoordinate* d_ac, int n_atoms,
                              char* d_lines, cudaStream_t stream);

// Format a single AtomCoordinate into a PDB_ATOM_LINE_LEN-byte line on the host
// (used to patch OXT slots and raw/CPU-fallback fragments byte-identically).
void formatAtomLineHost(const AtomCoordinate& atom, char* dst);

#endif // FOLDCOMP_WITH_CUDA
