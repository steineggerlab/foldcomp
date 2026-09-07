// SPDX-License-Identifier: MIT
/**
 * File: bond_constants.h
 * Project: foldcomp
 * Description:
 *     Idealized backbone bond-length constants (Angstroms) shared by the CPU
 *     NeRF reconstruction (foldcomp.h/.cpp) and the GPU NeRF kernel
 *     (gpu_nerf.cu), kept in one place so the two implementations can't drift.
 */
#pragma once

#define N_TO_CA_DIST 1.4581
#define CA_TO_C_DIST 1.5281
#define C_TO_N_DIST 1.3311
#define PRO_N_TO_CA_DIST 1.353
