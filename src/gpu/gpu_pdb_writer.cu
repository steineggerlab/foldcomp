// SPDX-License-Identifier: MIT
/**
 * File: gpu_pdb_writer.cu
 * Project: foldcomp
 * Description: GPU kernel + host launcher for formatting PDB ATOM records.
 */
#ifdef FOLDCOMP_WITH_CUDA

#include "gpu_pdb_writer.h"
#include "pdb_format.h"

#include <cuda_runtime.h>

// One thread per atom: writes a fixed PDB_ATOM_LINE_LEN-byte ATOM record.
__global__ void format_pdb_atom_lines_kernel(
    const AtomCoordinate* __restrict__ d_ac, int n_atoms, char* __restrict__ d_lines) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_atoms) return;
    pf_format_atom_line(d_ac[i], d_lines + static_cast<size_t>(i) * PDB_ATOM_LINE_LEN);
}

void formatAtomLinesDeviceGPU(const AtomCoordinate* d_ac, int n_atoms,
                              char* d_lines, cudaStream_t stream) {
    if (n_atoms <= 0) return;
    constexpr int BLOCK = 256;
    int grid = (n_atoms + BLOCK - 1) / BLOCK;
    format_pdb_atom_lines_kernel<<<grid, BLOCK, 0, stream>>>(
        d_ac, n_atoms, d_lines);
}

void formatAtomLineHost(const AtomCoordinate& atom, char* dst) {
    pf_format_atom_line(atom, dst);
}

#endif // FOLDCOMP_WITH_CUDA
