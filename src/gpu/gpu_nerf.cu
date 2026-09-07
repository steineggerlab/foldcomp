// SPDX-License-Identifier: MIT

/**
 * File: gpu_nerf.cu
 * Project: foldcomp
 * Created: 2026-01-06
 * Description:
 *     GPU-accelerated NeRF (Natural Extension of Reference Frame) operations using CUDA.
 *     Implements the place_atom algorithm for parallel atom coordinate calculation.
 */

#include "gpu_nerf.h"
#include "bond_constants.h"
#include "gpu_context.h"
#include "gpu_buffer.h"
#include "gpu_stream.h"
#include <nvtx3/nvtx3.hpp>
#include <cuda_runtime.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <cstdint>

static constexpr int8_t PRO_AA_IDX = 14; // AA_PRO = 14

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Device helper functions for vector operations
__device__ inline float3 cross_product_device(const float3& a, const float3& b) {
    return make_float3(
        a.y * b.z - b.y * a.z,
        a.z * b.x - b.z * a.x,
        a.x * b.y - b.x * a.y);
}

__device__ inline float norm_device(const float3& v) {
    return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}

/**
 * Place one atom via NeRF given three reference atoms (a,b,c) and bond params.
 * Angles are in radians. Returns the placed atom's coordinates.
 */
__device__ inline float3 nerf_place(
    const float3& a, const float3& b, const float3& c,
    float bond_length, float bond_angle, float torsion_angle) {
    float3 ab = make_float3(b.x - a.x, b.y - a.y, b.z - a.z);
    float3 bc = make_float3(c.x - b.x, c.y - b.y, c.z - b.z);

    float bc_norm = norm_device(bc);
    float3 bcn = make_float3(bc.x / bc_norm, bc.y / bc_norm, bc.z / bc_norm);

    float curr_x = -bond_length * cosf(bond_angle);
    float curr_y = bond_length * cosf(torsion_angle) * sinf(bond_angle);
    float curr_z = bond_length * sinf(torsion_angle) * sinf(bond_angle);

    float3 n = cross_product_device(ab, bcn);
    float n_norm = norm_device(n);
    n.x /= n_norm;
    n.y /= n_norm;
    n.z /= n_norm;

    float3 nbc = cross_product_device(n, bcn);

    return make_float3(
        bcn.x * curr_x + nbc.x * curr_y + n.x * curr_z + c.x,
        bcn.y * curr_x + nbc.y * curr_y + n.y * curr_z + c.y,
        bcn.z * curr_x + nbc.z * curr_y + n.z * curr_z + c.z);
}

// Geometric bond angle (degrees) at `b`, given its two neighbors — matches the
// CPU's float3d.h angle(atm1, atm2, atm3): atan2(|cross|, dot) of the two
// vectors from b, which is well-conditioned for angles near 0/180 (unlike
// acos(dot/|.||.|)).
__device__ inline float angle_device(const float3& a, const float3& b, const float3& c) {
    float3 d1 = make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
    float3 d2 = make_float3(c.x - b.x, c.y - b.y, c.z - b.z);
    float3 cr = cross_product_device(d1, d2);
    float cross_norm = norm_device(cr);
    float dot = d1.x * d2.x + d1.y * d2.y + d1.z * d2.z;
    return atan2f(cross_norm, dot) * 180.0f / (float)M_PI;
}

// Bond length constants: N_TO_CA_DIST/CA_TO_C_DIST/C_TO_N_DIST/PRO_N_TO_CA_DIST
// come from bond_constants.h (shared with foldcomp.h) above.

// ============================================================================
// Anchor-segment parallel backbone reconstruction
// ============================================================================

/**
 * Initialize per-unit (protein × segment) sliding-window prev_atoms buffers
 * from the protein-level prevAtoms and stored anchor coordinates.
 *
 * Unit index = protein_idx * max_n_segments + segment_idx
 *
 * Segment 0 of protein p uses prevAtoms_flat[p*9].
 * Segment s≥1 uses anchorCoordinates_flat[anchor_coord_offsets[p] + (s-1)*9].
 */
__global__ void init_segment_prev_atoms_kernel(
    float* d_unit_prev_atoms,          // [total_units * 9], output
    const float* d_init_prev_atoms,    // [n_proteins * 9]
    const float* d_anchor_coords,      // flat anchor coord data (floats)
    const int* d_anchor_coord_offsets, // [n_proteins+1], float indices into d_anchor_coords
    const int* d_n_anchors,            // [n_proteins]
    int n_proteins,
    int max_n_segments // stride = max(n_anchors[p]-1)
) {
    int unit_idx = blockIdx.x * blockDim.x + threadIdx.x;
    int protein_idx = unit_idx / max_n_segments;
    int segment_idx = unit_idx % max_n_segments;

    if (protein_idx >= n_proteins) return;
    if (segment_idx >= d_n_anchors[protein_idx] - 1) return;

    float* dst = d_unit_prev_atoms + unit_idx * 9;

    if (segment_idx == 0) {
        const float* src = d_init_prev_atoms + protein_idx * 9;
        for (int i = 0; i < 9; i++) dst[i] = src[i];
    } else {
        // anchor_coord_offsets is a float index; (s-1)*9 gives inner anchor s-1
        int float_off = d_anchor_coord_offsets[protein_idx] + (segment_idx - 1) * 9;
        const float* src = d_anchor_coords + float_off;
        for (int i = 0; i < 9; i++) dst[i] = src[i];
    }
}

/**
 * Fused segment kernel: one thread processes its entire anchor segment (all residues,
 * all 3 waves) in a single launch.  Eliminates 3×max_seg_len separate kernel launches
 * by keeping the sliding window in registers instead of round-tripping through
 * d_unit_prev_atoms every wave.
 *
 * Replaces the outer (local_res × wave) loop on the CPU side with an inner
 * per-thread loop.  Number of launches: 1 (vs. 3 × max_seg_len before).
 */
__global__ void backbone_segment_full_kernel(
    const float* d_unit_prev_atoms, // [total_units * 9], read-only (initial anchor state)
    const float* d_angles,
    const int* d_angle_offsets,
    const int8_t* d_residue_types,
    const int* d_residue_offsets,
    const int* d_anchor_indices,
    const int* d_anchor_offsets,
    const int* d_n_anchors,
    const int* d_n_residues,
    float* d_output,
    const int* d_coord_offsets,
    int n_proteins,
    int max_n_segments) {
    int unit_idx = blockIdx.x * blockDim.x + threadIdx.x;
    int protein_idx = unit_idx / max_n_segments;
    int segment_idx = unit_idx % max_n_segments;

    if (protein_idx >= n_proteins) return;

    int n_anchors = d_n_anchors[protein_idx];
    if (segment_idx >= n_anchors - 1) return;

    int anchor_base = d_anchor_offsets[protein_idx];
    int anchor_start = d_anchor_indices[anchor_base + segment_idx];
    int anchor_end = d_anchor_indices[anchor_base + segment_idx + 1];
    int seg_len = anchor_end - anchor_start;
    int n_res = d_n_residues[protein_idx];
    int res_base = d_residue_offsets[protein_idx];
    int angle_base_p = d_angle_offsets[protein_idx];
    int out_base_p = d_coord_offsets[protein_idx];

    // Load initial sliding-window state from anchor (read-only; no write-back needed)
    const float* my_prev = d_unit_prev_atoms + unit_idx * 9;
    float3 a = make_float3(my_prev[0], my_prev[1], my_prev[2]);
    float3 b = make_float3(my_prev[3], my_prev[4], my_prev[5]);
    float3 c = make_float3(my_prev[6], my_prev[7], my_prev[8]);

    const float k = (float)(M_PI / 180.0);

    for (int local_res = 0; local_res < seg_len; local_res++) {
        int global_res = anchor_start + local_res;
        if (global_res >= n_res - 1) break;

        int angle_base = angle_base_p + global_res * 6;
        float phi = d_angles[angle_base + 0];
        float psi = d_angles[angle_base + 1];
        float omega = d_angles[angle_base + 2];
        float n_ca_c_angle = d_angles[angle_base + 3];
        float ca_c_n_angle = d_angles[angle_base + 4];
        float c_n_ca_angle = d_angles[angle_base + 5];
        int8_t next_res = d_residue_types[res_base + global_res]; // current residue, matching CPU behaviour
        int out_off = out_base_p + global_res * 9;

        // Wave 0: N
        float3 d = nerf_place(a, b, c, C_TO_N_DIST, ca_c_n_angle * k, psi * k);
        d_output[out_off + 0] = d.x;
        d_output[out_off + 1] = d.y;
        d_output[out_off + 2] = d.z;
        a = b;
        b = c;
        c = d;

        // Wave 1: CA
        float bl_ca = (next_res == PRO_AA_IDX) ? PRO_N_TO_CA_DIST : N_TO_CA_DIST;
        d = nerf_place(a, b, c, bl_ca, c_n_ca_angle * k, omega * k);
        d_output[out_off + 3] = d.x;
        d_output[out_off + 4] = d.y;
        d_output[out_off + 5] = d.z;
        a = b;
        b = c;
        c = d;

        // Wave 2: C
        d = nerf_place(a, b, c, CA_TO_C_DIST, n_ca_c_angle * k, phi * k);
        d_output[out_off + 6] = d.x;
        d_output[out_off + 7] = d.y;
        d_output[out_off + 8] = d.z;
        a = b;
        b = c;
        c = d;
    }
    // Sliding window state stays in registers — no write-back to d_unit_prev_atoms needed.
}

// Forward-reconstructed position at local atom index m within a segment (m<3 is the
// segment's seed atom, from d_unit_prev_atoms; m>=3 was written by
// backbone_segment_full_kernel at global residue index anchor_start+(m-3)/3 — same
// indexing that kernel used: out_off = out_base_p + global_res*9).
__device__ inline float3 segment_fwd_atom(
    const float* d_unit_prev_atoms, const float* d_output,
    int unit_idx, int out_base_p, int anchor_start, int m) {
    if (m < 3) {
        const float* p = d_unit_prev_atoms + unit_idx * 9 + m * 3;
        return make_float3(p[0], p[1], p[2]);
    }
    int mm = m - 3;
    int off = out_base_p + (anchor_start + mm / 3) * 9 + (mm % 3) * 3;
    return make_float3(d_output[off], d_output[off + 1], d_output[off + 2]);
}

// Blend forward position f with backward position b (weighted by distance from each
// segment end, matching the CPU's (fwd*(total-m) + bwd*m)/total) and write it into
// d_output at local atom index m (m must be >= 3 — the seed atoms are never written).
__device__ inline void segment_write_blended(
    float* d_output, int out_base_p, int anchor_start, int total, int m,
    const float3& f, const float3& b) {
    const float wf = (float)(total - m), wb = (float)m;
    int mm = m - 3;
    int off = out_base_p + (anchor_start + mm / 3) * 9 + (mm % 3) * 3;
    d_output[off + 0] = (f.x * wf + b.x * wb) / total;
    d_output[off + 1] = (f.y * wf + b.y * wb) / total;
    d_output[off + 2] = (f.z * wf + b.z * wb) / total;
}

// ----------------------------------------------------------------------------
// Rigid-frame-transfer correction (pass 2)
//
// NeRF is exactly equivariant under rigid transforms: nerf_place(a,b,c,...) =
// c + R(a,b,c)*v, where R = [bcn | nbc | n] depends only on the *directions*
// (b-a) and (c-b), and v depends only on the bond parameters. So replacing a
// segment's seed triple (a,b,c) with a rigidly-transformed (Q*a+t, Q*b+t,
// Q*c+t) transforms every downstream atom the segment computes by that exact
// same (Q,t), by induction on the recursive sliding window.
//
// Frame3 below is that same [bcn | nbc | n] rotation (columns), with origin
// at c, matching nerf_place/angle_device's convention exactly.
// ----------------------------------------------------------------------------
struct Frame3 {
    float3 origin, bcn, nbc, n;
};

__device__ inline Frame3 build_frame_device(const float3& a, const float3& b, const float3& c) {
    float3 ab = make_float3(b.x - a.x, b.y - a.y, b.z - a.z);
    float3 bc = make_float3(c.x - b.x, c.y - b.y, c.z - b.z);
    float bc_norm = norm_device(bc);
    float3 bcn = make_float3(bc.x / bc_norm, bc.y / bc_norm, bc.z / bc_norm);

    float3 n = cross_product_device(ab, bcn);
    float n_norm = norm_device(n);
    n.x /= n_norm;
    n.y /= n_norm;
    n.z /= n_norm;

    float3 nbc = cross_product_device(n, bcn);

    Frame3 f;
    f.origin = c;
    f.bcn = bcn;
    f.nbc = nbc;
    f.n = n;
    return f;
}

// Express world-space point p in frame f's local (orthonormal) coordinates.
__device__ inline float3 frame_to_local(const Frame3& f, const float3& p) {
    float3 d = make_float3(p.x - f.origin.x, p.y - f.origin.y, p.z - f.origin.z);
    return make_float3(
        d.x * f.bcn.x + d.y * f.bcn.y + d.z * f.bcn.z,
        d.x * f.nbc.x + d.y * f.nbc.y + d.z * f.nbc.z,
        d.x * f.n.x + d.y * f.n.y + d.z * f.n.z);
}

// Reconstruct a world-space point from frame f's local coordinates (inverse of
// frame_to_local, valid since bcn/nbc/n are orthonormal).
__device__ inline float3 frame_from_local(const Frame3& f, const float3& loc) {
    return make_float3(
        f.bcn.x * loc.x + f.nbc.x * loc.y + f.n.x * loc.z + f.origin.x,
        f.bcn.y * loc.x + f.nbc.y * loc.y + f.n.y * loc.z + f.origin.y,
        f.bcn.z * loc.x + f.nbc.z * loc.y + f.n.z * loc.z + f.origin.z);
}

/**
 * Per-segment backward NeRF pass + forward/backward blend — the GPU counterpart of
 * the CPU's reconstructBackboneReverse/Nerf::reconstructWithReversed (foldcomp.cpp).
 * Forward-only reconstruction lets error accumulate across a whole segment; the CPU
 * bounds it by also walking backward from the segment's known-exact end anchor and
 * blending (fwd*(total-i) + bwd*i)/total, so both ends are pinned and error peaks at
 * the segment midpoint instead of growing unbounded toward the end. Must run after
 * backbone_segment_full_kernel (same stream) since it reads that kernel's output.
 *
 * Runs entirely from data this segment already owns (its own forward output in
 * d_output, its own seed in d_unit_prev_atoms, its own end anchor in
 * d_anchor_coords) — no cross-segment dependency, so it stays embarrassingly
 * parallel like the forward kernel.
 *
 * Local atom numbering within a segment (matching the CPU's `atom` array in
 * reconstructBackboneReverse): m=0..2 is the segment's seed residue (already exact
 * or CPU-equivalent — untouched here), m=3..total-1 is the seg_len freshly-computed
 * residues, where total=(seg_len+1)*3. m=total-3..total-1 is the boundary residue,
 * whose true value comes from d_anchor_coords, not NeRF.
 */
__global__ void backbone_segment_backward_blend_kernel(
    float* d_output,
    const float* d_unit_prev_atoms,
    const float* d_angles,
    const int* d_angle_offsets,
    const int* d_anchor_indices,
    const int* d_anchor_offsets,
    const float* d_anchor_coords,
    const int* d_anchor_coord_offsets,
    const int* d_n_anchors,
    const int* d_coord_offsets,
    int n_proteins,
    int max_n_segments) {
    int unit_idx = blockIdx.x * blockDim.x + threadIdx.x;
    int protein_idx = unit_idx / max_n_segments;
    int segment_idx = unit_idx % max_n_segments;
    if (protein_idx >= n_proteins) return;

    int n_anchors = d_n_anchors[protein_idx];
    if (segment_idx >= n_anchors - 1) return;

    int anchor_base = d_anchor_offsets[protein_idx];
    int anchor_start = d_anchor_indices[anchor_base + segment_idx];
    int anchor_end = d_anchor_indices[anchor_base + segment_idx + 1];
    int seg_len = anchor_end - anchor_start;
    if (seg_len < 1) return; // degenerate segment: nothing to blend

    int out_base_p = d_coord_offsets[protein_idx];
    int angle_base_p = d_angle_offsets[protein_idx];
    const int total = (seg_len + 1) * 3;

    // This segment's true end-anchor coordinates (N, CA, C of the boundary residue) —
    // stored losslessly at compression time, unlike everything reconstructed via NeRF.
    int aoff = d_anchor_coord_offsets[protein_idx] + segment_idx * 9;
    const float3 exact_N  = make_float3(d_anchor_coords[aoff + 0], d_anchor_coords[aoff + 1], d_anchor_coords[aoff + 2]);
    const float3 exact_CA = make_float3(d_anchor_coords[aoff + 3], d_anchor_coords[aoff + 4], d_anchor_coords[aoff + 5]);
    const float3 exact_C  = make_float3(d_anchor_coords[aoff + 6], d_anchor_coords[aoff + 7], d_anchor_coords[aoff + 8]);

    // Seed the backward sliding window from the boundary residue's exact coords —
    // bwd[total-1]=C, bwd[total-2]=CA, bwd[total-3]=N — and blend+write those three
    // positions immediately (their "backward" value needs no NeRF step). Capture the
    // pure forward reads before write_blended overwrites d_output at these positions:
    // angle_device below needs the pristine forward geometry (matching CPU's
    // getBondAngles running on a `forwardAtom` snapshot taken before any blending),
    // not the already-blended value that a later fwd(m) read would return.
    const float3 f_tm3 = segment_fwd_atom(d_unit_prev_atoms, d_output, unit_idx, out_base_p, anchor_start, total - 3);
    const float3 f_tm2 = segment_fwd_atom(d_unit_prev_atoms, d_output, unit_idx, out_base_p, anchor_start, total - 2);
    const float3 f_tm1 = segment_fwd_atom(d_unit_prev_atoms, d_output, unit_idx, out_base_p, anchor_start, total - 1);
    segment_write_blended(d_output, out_base_p, anchor_start, total, total - 3, f_tm3, exact_N);
    segment_write_blended(d_output, out_base_p, anchor_start, total, total - 2, f_tm2, exact_CA);
    segment_write_blended(d_output, out_base_p, anchor_start, total, total - 1, f_tm1, exact_C);
    float3 p1 = exact_N, p2 = exact_CA, p3 = exact_C; // bwd[m+1], bwd[m+2], bwd[m+3] for m=total-4
    float3 pf1 = f_tm3, pf2 = f_tm2;                  // pure fwd[m+1], fwd[m+2] for m=total-4

    const float k_deg = (float)(M_PI / 180.0);

    for (int m = total - 4; m >= 3; m--) {
        const float3 pf0 = segment_fwd_atom(d_unit_prev_atoms, d_output, unit_idx, out_base_p, anchor_start, m); // pristine forward read — not yet blended
        const int typ = m % 3; // 0=N, 1=CA, 2=C (position within its residue)
        // gres = the torsion/bond-angle *record* governing this atom's local residue
        // K=m/3: forward uses record (K-1) to *place* residue K from residue (K-1) —
        // that's segment_fwd_atom's read location — but the CPU's backward pass reuses
        // the *outgoing* record K (the one that places residue K+1 from residue K) for
        // all 3 atoms of residue K instead (verified empirically against the CPU's own
        // reconstructWithReversed — this asymmetry is inherent to that algorithm, not a
        // derivable "should"). record index is a whole-protein index, not segment-local.
        const int gres = anchor_start + m / 3;

        float bond_length;
        if (typ == 2) {
            bond_length = C_TO_N_DIST;
        } else if (typ == 1) {
            bond_length = CA_TO_C_DIST;
        } else {
            bond_length = N_TO_CA_DIST;
        }

        // Geometric bond angle at m+1, from the forward reconstruction's own
        // geometry (not the discretized/continuized value) — matches CPU's
        // getBondAngles(forwardAtom) fed into reconstructWithReversed.
        const float bond_angle = angle_device(pf0, pf1, pf2);

        const int abase = angle_base_p + gres * 6;
        // torsion(m): N uses psi, CA uses omega, C uses phi — same assignment
        // backbone_segment_full_kernel uses to forward-place each atom.
        const float torsion = (typ == 0) ? d_angles[abase + 1]
                             : (typ == 1) ? d_angles[abase + 2]
                                          : d_angles[abase + 0];

        float3 bwd_m = nerf_place(p3, p2, p1, bond_length, bond_angle * k_deg, torsion * k_deg);
        segment_write_blended(d_output, out_base_p, anchor_start, total, m, pf0, bwd_m);

        p3 = p2;
        p2 = p1;
        p1 = bwd_m;
        pf2 = pf1;
        pf1 = pf0;
    }
}

/**
 * Rigid-frame-transfer correction pass — cheaply approximates what the CPU actually
 * does (seed every segment ≥1 from the *previous* segment's own post-blend boundary
 * result, a serial dependency — see prevForAnchor in Foldcomp::decompressBackbone,
 * foldcomp.cpp ~line 1015-1053) without re-running any NeRF/trig work.
 *
 * backbone_segment_full_kernel/backbone_segment_backward_blend_kernel above already
 * ran this segment's forward+backward+blend using the raw *exact* anchor coordinate
 * as its seed (loaded once, up front, by init_segment_prev_atoms_kernel into
 * d_unit_prev_atoms — embarrassingly parallel, no cross-segment dependency). That
 * seed is slightly "wrong" vs. what the CPU would use: the CPU's own reconstruction
 * of the same boundary residue is instead the *blended* forward/backward result from
 * the previous segment, which differs very slightly from the raw exact anchor.
 *
 * Since NeRF is exactly equivariant under rigid transforms (nerf_place(a,b,c,...) =
 * c + R(a,b,c)*v — see Frame3 doc above), we don't need to redo any trig: we just
 * need the single rigid transform (Q,t) that carries the OLD seed frame (raw exact
 * anchor, still sitting untouched in d_unit_prev_atoms) onto the NEW seed frame (the
 * previous segment's already-blended boundary result, now sitting in d_output — see
 * indexing note below), and apply that same transform to every atom this segment
 * already computed. This is an O(seg_len) loop of dot products/linear combinations —
 * no re-blending, no torsion/bond-angle lookups, no re-running kernels 2/3.
 *
 * The two seed triples (raw exact anchor vs. CPU-equivalent blended value) are not
 * exactly congruent — blending introduces tiny numerical differences — so the rigid
 * transform is a best fit in spirit but computed exactly from the 3 corresponding
 * points; since both triples describe near-identical bond geometry for the same
 * physical residue, this is an excellent approximation of what a full CPU-accurate
 * reseed+rerun would produce, at a fraction of the cost.
 *
 * Must run after backbone_segment_backward_blend_kernel (same stream): it reads the
 * *previous* segment's finished, blended d_output slot as this segment's new seed.
 *
 * Indexing: this segment's own seed residue is global index anchor_start (the exact
 * anchor at which this segment begins — see init_segment_prev_atoms_kernel's use of
 * anchor_coords[(segment_idx-1)*9]). backbone_segment_full_kernel/backward_blend_kernel
 * write global residue g into d_output slot out_base_p+g*9 for slot content = residue
 * g+1 (out_off=out_base_p+global_res*9, global_res ranges anchor_start..anchor_end-1,
 * writing residues anchor_start+1..anchor_end). The *previous* segment's last written
 * slot has global_res=anchor_start-1 (since its anchor_end == this segment's
 * anchor_start), i.e. content = residue anchor_start — exactly this segment's own
 * seed residue, post-blend. So the NEW seed lives at out_base_p+(anchor_start-1)*9.
 *
 * Segment 0 is skipped: its seed is the protein's own true starting atoms (CPU-exact
 * already), not an approximation to correct.
 *
 * Correction strength fades toward the segment's far end (see the m-dependent weight
 * w = (total-m)/total in the write loop below) — this weight is *derived*, not fitted.
 * backbone_segment_backward_blend_kernel's blended output at local atom m is
 * p(m) = (fwd(m)*(total-m) + bwd(m)*m)/total, where fwd is seed-dependent (the quantity
 * a seed correction Q needs to fix) and bwd is seed-*independent* (walked backward from
 * this segment's own exact end anchor, no dependency on the start seed at all).
 * Correcting only the seed gives the exact corrected target:
 *   p'(m) = (Q·fwd(m)*(total-m) + bwd(m)*m)/total
 *         = p(m) + (total-m)/total * (Q·fwd(m) - fwd(m))
 * i.e. exactly w = (total-m)/total applied to the *forward* displacement
 * Q·fwd(m) - fwd(m). Applying the transform at full strength (unweighted, w=1
 * everywhere) isn't an empirically-discovered mistake — it applies a seed correction to
 * the m/total share of p(m) that has no seed dependence at all, which is exactly why it
 * measurably increases RMSD instead of reducing it.
 *
 * Residual approximation: the write loop below applies Q to the already-blended p(m)
 * (read back from d_output) rather than to the pure fwd(m) term, i.e. it computes
 * w*(Q·p - p) instead of the exact w*(Q·fwd - fwd). Since Q is a near-identity rigid
 * transform, the gap between these is (R-I)*(m/total)*(bwd-fwd) — a near-identity
 * rotation applied to a sub-0.1 Å difference, order 1e-4 Å, well below the ~0.0064 Å
 * measured residual (see README_GPU.md's accuracy-divergence section) — so it is left
 * as the one remaining approximation rather than plumbing a separate fwd(m) value
 * through just to close it.
 */
__global__ void backbone_segment_frame_transfer_kernel(
    float* d_output,
    const float* d_unit_prev_atoms,
    const int* d_anchor_indices,
    const int* d_anchor_offsets,
    const int* d_n_anchors,
    const int* d_coord_offsets,
    int n_proteins,
    int max_n_segments) {
    int unit_idx = blockIdx.x * blockDim.x + threadIdx.x;
    int protein_idx = unit_idx / max_n_segments;
    int segment_idx = unit_idx % max_n_segments;
    if (protein_idx >= n_proteins) return;

    int n_anchors = d_n_anchors[protein_idx];
    if (segment_idx >= n_anchors - 1) return;
    if (segment_idx == 0) return; // seed already CPU-exact; nothing to correct

    int anchor_base = d_anchor_offsets[protein_idx];
    int anchor_start = d_anchor_indices[anchor_base + segment_idx];
    int anchor_end = d_anchor_indices[anchor_base + segment_idx + 1];
    int seg_len = anchor_end - anchor_start;
    if (seg_len < 1) return; // degenerate segment: nothing was written by pass 1

    int out_base_p = d_coord_offsets[protein_idx];
    const int total = (seg_len + 1) * 3;

    // OLD frame: raw exact-anchor seed this segment's pass-1 kernels actually used
    // (untouched by any of the earlier kernels).
    const float* old_p = d_unit_prev_atoms + unit_idx * 9;
    float3 oa = make_float3(old_p[0], old_p[1], old_p[2]);
    float3 ob = make_float3(old_p[3], old_p[4], old_p[5]);
    float3 oc = make_float3(old_p[6], old_p[7], old_p[8]);

    // NEW frame: the previous segment's post-blend result for this same physical
    // boundary residue (see indexing note above).

    if (anchor_start < 1) return; //  no previous segment output to reseed from
    int new_off = out_base_p + (anchor_start - 1) * 9;
    float3 na = make_float3(d_output[new_off + 0], d_output[new_off + 1], d_output[new_off + 2]);
    float3 nb = make_float3(d_output[new_off + 3], d_output[new_off + 4], d_output[new_off + 5]);
    float3 nc = make_float3(d_output[new_off + 6], d_output[new_off + 7], d_output[new_off + 8]);

    Frame3 old_frame = build_frame_device(oa, ob, oc);
    Frame3 new_frame = build_frame_device(na, nb, nc);

    // Re-express this segment's own *interior* atoms (m=3..total-4) in OLD-frame local
    // coordinates, then reconstruct them in NEW-frame world coordinates. Deliberately
    // NOT m=3..total-1 (the full range backbone_segment_backward_blend_kernel wrote):
    // m=total-3..total-1 is this segment's own END-boundary residue (global index
    // anchor_end), which is 1) already pinned near-exact to this segment's *own* true
    // end anchor by the blend kernel — unrelated to the START-boundary mismatch this
    // transform corrects, so rotating it by this Q would only pull it away from that
    // exact value — and 2) is the *next* segment's NEW seed (read from this same
    // d_output slot, computed above via the identical (anchor_start-1) indexing
    // relation), which this same kernel launch processes in parallel with no ordering
    // guarantee — touching it here would race with that read. So it must stay exactly
    // as backbone_segment_backward_blend_kernel left it, both for correctness and to
    // avoid a cross-thread hazard.
    for (int m = 3; m < total - 3; m++) {
        int mm = m - 3;
        int off = out_base_p + (anchor_start + mm / 3) * 9 + (mm % 3) * 3;
        float3 p = make_float3(d_output[off + 0], d_output[off + 1], d_output[off + 2]);
        float3 loc = frame_to_local(old_frame, p);
        float3 q = frame_from_local(new_frame, loc); // fully-corrected (seed-transformed) position

        // Blend p->q by the same (total-m)/total weight backbone_segment_backward_blend_kernel
        // used for its own forward/backward blend — see kernel doc above for why.
        float w = (float)(total - m) / (float)total;
        d_output[off + 0] = p.x + w * (q.x - p.x);
        d_output[off + 1] = p.y + w * (q.y - p.y);
        d_output[off + 2] = p.z + w * (q.z - p.z);
    }
}

void reconstructBackboneGPU_segmented(
    NerfBuffers& ctx,
    int n_proteins,
    const std::vector<int>& n_residues_per_protein,
    const std::vector<float>& prev_atoms_flat, // [n_proteins * 9]
    const std::vector<int>& backbone_offsets,  // [n_proteins+1] residue counts (angle offsets = *6)
    const std::vector<int8_t>& residue_types,
    const std::vector<int>& residue_offsets,       // [n_proteins+1]
    const std::vector<int>& anchor_indices_flat,   // batch.anchorIndices_flat
    const std::vector<int>& anchor_offsets,        // batch.anchor_offsets [n_proteins+1]
    const std::vector<float>& anchor_coords_flat,  // batch.anchorCoordinates_flat
    const std::vector<int>& anchor_coord_offsets,  // batch.anchor_coord_offsets [n_proteins+1]
    const std::vector<int>& n_anchors_per_protein, // batch.nAllAnchor_batch
    PinnedFloatVec& output_coords,
    std::vector<int>& coord_offsets,
    const float* d_angles,
    cudaStream_t stream) {
    if (n_proteins == 0) {
        output_coords.clear();
        coord_offsets.clear();
        return;
    }

    coord_offsets.resize(n_proteins + 1);
    coord_offsets[0] = 0;
    for (int i = 0; i < n_proteins; i++) {
        coord_offsets[i + 1] = coord_offsets[i] + (n_residues_per_protein[i] - 1) * 9;
    }
    int total_output = coord_offsets[n_proteins];

    // --- Derive loop bounds ---
    int max_n_anchors, max_n_segments, total_units, max_seg_len;
    max_n_anchors = *std::max_element(n_anchors_per_protein.begin(), n_anchors_per_protein.end());
    max_n_segments = max_n_anchors - 1;
    total_units = n_proteins * max_n_segments;

    if (total_units <= 0 && total_output > 0) {
        throw std::runtime_error(
            "reconstructBackboneGPU_segmented: no valid anchor segments (max_n_anchors=" +
            std::to_string(max_n_anchors) + ") but " + std::to_string(total_output) +
            " output floats expected — malformed anchor data");
    }

    max_seg_len = 0;
    for (int p = 0; p < n_proteins; p++) {
        int ab = anchor_offsets[p];
        for (int s = 0; s < n_anchors_per_protein[p] - 1; s++) {
            int gap = anchor_indices_flat[ab + s + 1] - anchor_indices_flat[ab + s];
            if (gap > max_seg_len) max_seg_len = gap;
        }
    }

    // Compute pack offsets (all 4-byte unless noted; rcodes is int8_t at the end)
    const size_t sz_prev = (size_t)n_proteins * 9 * sizeof(float);
    const size_t sz_aoff = backbone_offsets.size() * sizeof(int); // angle offsets derived as *6
    const size_t sz_roff = residue_offsets.size() * sizeof(int);
    const size_t sz_aidx = anchor_indices_flat.size() * sizeof(int);
    const size_t sz_aoffg = anchor_offsets.size() * sizeof(int);
    const size_t sz_acoords = anchor_coords_flat.size() * sizeof(float);
    const size_t sz_acoff = anchor_coord_offsets.size() * sizeof(int);
    const size_t sz_nanch = (size_t)n_proteins * sizeof(int);
    const size_t sz_nres = (size_t)n_proteins * sizeof(int);
    const size_t sz_coff = coord_offsets.size() * sizeof(int);
    const size_t sz_rcodes = residue_types.size() * sizeof(int8_t);

    const size_t off_prev = 0;
    const size_t off_aoff = off_prev + sz_prev;
    const size_t off_roff = off_aoff + sz_aoff;
    const size_t off_aidx = off_roff + sz_roff;
    const size_t off_aoffg = off_aidx + sz_aidx;
    const size_t off_acoords = off_aoffg + sz_aoffg;
    const size_t off_acoff = off_acoords + sz_acoords;
    const size_t off_nanch = off_acoff + sz_acoff;
    const size_t off_nres = off_nanch + sz_nanch;
    const size_t off_coff = off_nres + sz_nres;
    const size_t off_rcodes = off_coff + sz_coff; // int8_t: no alignment gap needed at end
    const size_t pack_sz = off_rcodes + sz_rcodes;

    char* h_pack = static_cast<char*>(ctx.ph_seg_pack.ensure(pack_sz));
    char* d_pack = static_cast<char*>(ctx.seg_pack.ensure(pack_sz, stream));

    float* h_prev = reinterpret_cast<float*>(h_pack + off_prev);
    int* h_aoff = reinterpret_cast<int*>(h_pack + off_aoff);
    int* h_roff = reinterpret_cast<int*>(h_pack + off_roff);
    int* h_aidx = reinterpret_cast<int*>(h_pack + off_aidx);
    int* h_aoffg = reinterpret_cast<int*>(h_pack + off_aoffg);
    float* h_acoords = reinterpret_cast<float*>(h_pack + off_acoords);
    int* h_acoff = reinterpret_cast<int*>(h_pack + off_acoff);
    int* h_nanch = reinterpret_cast<int*>(h_pack + off_nanch);
    int* h_nres = reinterpret_cast<int*>(h_pack + off_nres);
    int* h_coff = reinterpret_cast<int*>(h_pack + off_coff);
    int8_t* h_rcodes = reinterpret_cast<int8_t*>(h_pack + off_rcodes);

    float* d_init_prev_atoms = reinterpret_cast<float*>(d_pack + off_prev);
    int* d_angle_offsets = reinterpret_cast<int*>(d_pack + off_aoff);
    int* d_residue_offsets = reinterpret_cast<int*>(d_pack + off_roff);
    int* d_anchor_indices = reinterpret_cast<int*>(d_pack + off_aidx);
    int* d_anchor_offsets_gpu = reinterpret_cast<int*>(d_pack + off_aoffg);
    float* d_anchor_coords = reinterpret_cast<float*>(d_pack + off_acoords);
    int* d_anchor_coord_offsets_gpu = reinterpret_cast<int*>(d_pack + off_acoff);
    int* d_n_anchors = reinterpret_cast<int*>(d_pack + off_nanch);
    int* d_n_residues = reinterpret_cast<int*>(d_pack + off_nres);
    int* d_coord_offsets = reinterpret_cast<int*>(d_pack + off_coff);
    int8_t* d_residue_codes = reinterpret_cast<int8_t*>(d_pack + off_rcodes);
    ctx.last_d_rt = d_residue_codes;

    float* d_unit_prev_atoms = ctx.unit_prev_atoms.ensure_n<float>(total_units * 9, stream);
    float* d_output = ctx.output.ensure_n<float>(total_output, stream);

    memcpy(h_prev, prev_atoms_flat.data(), sz_prev);
    for (int i = 0; i < (int)backbone_offsets.size(); i++) h_aoff[i] = backbone_offsets[i] * 6;
    memcpy(h_rcodes, residue_types.data(), sz_rcodes);
    memcpy(h_roff, residue_offsets.data(), sz_roff);
    memcpy(h_aidx, anchor_indices_flat.data(), sz_aidx);
    memcpy(h_aoffg, anchor_offsets.data(), sz_aoffg);
    memcpy(h_acoords, anchor_coords_flat.data(), sz_acoords);
    memcpy(h_acoff, anchor_coord_offsets.data(), sz_acoff);
    memcpy(h_nanch, n_anchors_per_protein.data(), sz_nanch);
    memcpy(h_nres, n_residues_per_protein.data(), sz_nres);
    memcpy(h_coff, coord_offsets.data(), sz_coff);

    // --- Upload: single transfer for all per-protein metadata ---
    CUDA_CHECK(cudaMemcpyAsync(d_pack, h_pack, pack_sz, cudaMemcpyHostToDevice, stream));

    // --- Initialize per-unit prev_atoms from anchor data ---
    const int block_size = 256;
    int init_blocks = (total_units + block_size - 1) / block_size;
    init_segment_prev_atoms_kernel<<<init_blocks, block_size, 0, stream>>>(
        d_unit_prev_atoms,
        d_init_prev_atoms,
        d_anchor_coords,
        d_anchor_coord_offsets_gpu,
        d_n_anchors,
        n_proteins,
        max_n_segments);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
        throw std::runtime_error("init_segment_prev_atoms_kernel failed: " + std::string(cudaGetErrorString(err)));

    // Single fused launch: each thread processes its entire segment (all local_res × 3 waves)
    // in registers. Replaces 3 × max_seg_len separate kernel launches with 1.
    int seg_blocks = (total_units + block_size - 1) / block_size;
    backbone_segment_full_kernel<<<seg_blocks, block_size, 0, stream>>>(
        d_unit_prev_atoms,
        d_angles,
        d_angle_offsets,
        d_residue_codes,
        d_residue_offsets,
        d_anchor_indices,
        d_anchor_offsets_gpu,
        d_n_anchors,
        d_n_residues,
        d_output,
        d_coord_offsets,
        n_proteins,
        max_n_segments);

    // Backward pass + forward/backward blend: bounds intra-segment drift the same
    // way the CPU's reconstructBackboneReverse does (see kernel doc above). Must
    // follow the forward launch on the same stream since it reads d_output.
    backbone_segment_backward_blend_kernel<<<seg_blocks, block_size, 0, stream>>>(
        d_output,
        d_unit_prev_atoms,
        d_angles,
        d_angle_offsets,
        d_anchor_indices,
        d_anchor_offsets_gpu,
        d_anchor_coords,
        d_anchor_coord_offsets_gpu,
        d_n_anchors,
        d_coord_offsets,
        n_proteins,
        max_n_segments);
    err = cudaGetLastError();
    if (err != cudaSuccess)
        throw std::runtime_error(
            "backbone_segment_backward_blend_kernel failed: " + std::string(cudaGetErrorString(err)));

    // Rigid-frame-transfer correction: cheaply approximates the CPU's serial
    // segment-to-segment reseeding (see kernel doc above) by applying, to each
    // segment's already-computed atoms, the single rigid transform that carries its
    // raw-exact-anchor seed onto the previous segment's post-blend boundary result.
    // Must follow the backward/blend launch on the same stream (reads its output).
    backbone_segment_frame_transfer_kernel<<<seg_blocks, block_size, 0, stream>>>(
        d_output,
        d_unit_prev_atoms,
        d_anchor_indices,
        d_anchor_offsets_gpu,
        d_n_anchors,
        d_coord_offsets,
        n_proteins,
        max_n_segments);
    err = cudaGetLastError();
    if (err != cudaSuccess)
        throw std::runtime_error(
            "backbone_segment_frame_transfer_kernel failed: " + std::string(cudaGetErrorString(err)));

    // Backbone stays on device; expose pointer for sidechain kernels
    ctx.last_bb_dev_ptr = d_output;
}

// ============================================================================
// buildBackboneForSidechains: upload backbone metadata for sidechain kernels
// ============================================================================

BBSidechainPtrs buildBackboneForSidechains(
    NerfBuffers& ctx,
    int n_structs,
    const std::vector<int>& coord_offsets_cpu,
    const std::vector<int>& n_residues_per_protein,
    const std::vector<float>& prev_atoms_flat,
    int /*total_residues*/,
    cudaStream_t stream) {
    if (!ctx.last_bb_dev_ptr) return {nullptr, nullptr, nullptr, nullptr, nullptr};

    // Pack [float[n*9] pa | int[n+1] rfo | int[n+1] bbo] into a single H→D transfer.
    // Sidechain kernels use these pointers directly for inline backbone lookup.
    const size_t off_pa = 0;
    const size_t sz_pa = (size_t)n_structs * 9 * sizeof(float);
    const size_t off_rfo = sz_pa; // float and int both 4-byte — no gap
    const size_t sz_rfo = (size_t)(n_structs + 1) * sizeof(int);
    const size_t off_bbo = off_rfo + sz_rfo;
    const size_t sz_bbo = (size_t)(n_structs + 1) * sizeof(int);
    const size_t pack_sz = off_bbo + sz_bbo;

    char* h_pack = static_cast<char*>(ctx.ph_bb_sc_pack.ensure(pack_sz));
    char* d_pack = static_cast<char*>(ctx.d_bb_sc_pack.ensure(pack_sz, stream));

    float* h_pa = reinterpret_cast<float*>(h_pack + off_pa);
    int* h_rfo = reinterpret_cast<int*>(h_pack + off_rfo);
    int* h_bbo = reinterpret_cast<int*>(h_pack + off_bbo);

    h_rfo[0] = 0;
    for (int i = 0; i < n_structs; i++)
        h_rfo[i + 1] = h_rfo[i] + n_residues_per_protein[i];

    memcpy(h_pa, prev_atoms_flat.data(), sz_pa);
    memcpy(h_bbo, coord_offsets_cpu.data(), sz_bbo);

    CUDA_CHECK(cudaMemcpyAsync(d_pack, h_pack, pack_sz, cudaMemcpyHostToDevice, stream));

    return BBSidechainPtrs{
        ctx.last_bb_dev_ptr,
        reinterpret_cast<float*>(d_pack + off_pa),
        reinterpret_cast<int*>(d_pack + off_rfo),
        reinterpret_cast<int*>(d_pack + off_bbo),
        ctx.last_d_rt};
}
