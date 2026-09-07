// SPDX-License-Identifier: MIT
/**
 * File: gpu_foldcomp_batch.cpp
 * Project: foldcomp
 * Description:
 *     GPU batched-decompression implementations for BatchedDecompressionData
 *     and Foldcomp. Extracted from foldcomp.cpp.
 */
#include "foldcomp.h"

#include <stdexcept>

#include "amino_acid.h"
#include "utility.h"

#include <cmath>
#include <algorithm>
#include <cstdio>
#include <utility>

#include <nvtx3/nvtx3.hpp>
#include "gpu_discretizer.h"
#include "gpu_nerf.h"
#include "gpu_pdb_writer.h"
#include "gpu_sidechain.h"

void BatchedDecompressionData::finalize() {
    if (blobs.empty()) return;

    const size_t n = blobs.size();

    // Pre-compute totals for single-pass reservation
    size_t total_bb = 0, total_sc = 0, total_tf = 0;
    size_t total_ai = 0, total_ac = 0, total_res = 0;
    for (const auto& b : blobs) {
        total_bb  += (size_t)b.nResidue;
        total_sc  += (size_t)b.nSideChainTorsion;
        total_tf  += (size_t)b.nResidue;
        total_ai  += (size_t)b.nAllAnchor;
        total_ac  += (size_t)b.nAnchorGroups * 9;
        total_res += (size_t)b.nResidue;
    }

    backbone_offsets.reserve(n + 1);
    sidechain_offsets.reserve(n + 1);
    anchor_offsets.reserve(n + 1);
    anchor_coord_offsets.reserve(n + 1);
    residue_type_offsets.reserve(n + 1);
    nResidue_batch.reserve(n);
    nAtom_batch.reserve(n);
    nSideChainTorsion_batch.reserve(n);
    nAllAnchor_batch.reserve(n);
    firstResidue_batch.reserve(n);
    lastResidue_batch.reserve(n);
    hasOXT_batch.reserve(n);
    useAltAtomOrder_batch.reserve(n);
    titles_batch.reserve(n);
    tempFactorDisc_min_batch.reserve(n);
    tempFactorDisc_cont_batch.reserve(n);
    headers_batch.reserve(n);
    compressedBackBone_flat.reserve(total_bb);
    sideChainAnglesDiscretized_flat.reserve(total_sc);
    tempFactorsDiscretized_flat.reserve(total_tf);
    anchorIndices_flat.reserve(total_ai);
    anchorCoordinates_flat.reserve(total_ac);
    prevAtoms_flat.reserve(n * 9);
    OXT_coords_flat.reserve(n * 3);
    residues_flat.reserve(total_res);
    residue_types_flat.reserve(total_res);
    atom_offsets_flat.reserve(total_res + 1);  // +1 for the initial sentinel 0

    for (const auto& b : blobs) {
        // Per-structure metadata
        nResidue_batch.push_back(b.nResidue);
        nAtom_batch.push_back(b.nAtom);
        nSideChainTorsion_batch.push_back(b.nSideChainTorsion);
        nAllAnchor_batch.push_back(b.nAllAnchor);
        firstResidue_batch.push_back((char)b.firstResidue);
        lastResidue_batch.push_back((char)b.lastResidue);
        hasOXT_batch.push_back((char)b.hasOXT);
        useAltAtomOrder_batch.push_back(b.useAltAtomOrder);
        titles_batch.push_back(b.title);
        tempFactorDisc_min_batch.push_back(b.tempFactorDisc_min);
        tempFactorDisc_cont_batch.push_back(b.tempFactorDisc_cont);
        headers_batch.push_back(b.header);

        // Flat backbone raw bytes (pinned memory vector, 8 bytes/residue, file format)
        const uint64_t* bb = b.backbone_raw();
        compressedBackBone_flat.insert(compressedBackBone_flat.end(), bb, bb + b.nResidue);
        backbone_offsets.push_back((int)compressedBackBone_flat.size());

        // Flat sidechain angles (pinned memory vector)
        const uint8_t* sc = b.sc_angles();
        sideChainAnglesDiscretized_flat.insert(sideChainAnglesDiscretized_flat.end(),
                                               sc, sc + b.nSideChainTorsion);
        sidechain_offsets.push_back((int)sideChainAnglesDiscretized_flat.size());

        // Flat temp factors
        const uint8_t* tf = b.temp_factors();
        tempFactorsDiscretized_flat.insert(tempFactorsDiscretized_flat.end(),
                                           tf, tf + b.nResidue);

        // Flat anchor indices
        const int32_t* ai = b.anchor_idx();
        anchorIndices_flat.insert(anchorIndices_flat.end(), ai, ai + b.nAllAnchor);
        anchor_offsets.push_back((int)anchorIndices_flat.size());

        // Flat anchor coordinates
        const float* ac = b.anchor_coords();
        anchorCoordinates_flat.insert(anchorCoordinates_flat.end(),
                                      ac, ac + (size_t)b.nAnchorGroups * 9);
        anchor_coord_offsets.push_back((int)anchorCoordinates_flat.size());

        // Flat prevAtoms (9 floats per structure)
        const float* pa = b.prev_atoms();
        prevAtoms_flat.insert(prevAtoms_flat.end(), pa, pa + 9);

        // Flat OXT coords (3 floats per structure)
        const float* oxt = b.oxt_coords();
        OXT_coords_flat.insert(OXT_coords_flat.end(), oxt, oxt + 3);

        // Flat residues (one-letter char)
        const char* res = b.residues();
        residues_flat.insert(residues_flat.end(), res, res + b.nResidue);

        // Flat residue types (int8_t, AminoAcidIndex)
        const int8_t* rt = b.res_types();
        residue_types_flat.insert(residue_types_flat.end(), rt, rt + b.nResidue);
        residue_type_offsets.push_back((int)residue_types_flat.size());

        // Flat atom offsets: local (0-based) → global, one entry per residue
        const int32_t* ao = b.atom_offsets();
        for (int i = 0; i < b.nResidue; i++)
            atom_offsets_flat.push_back(_atom_running + ao[i]);
        _atom_running += ao[b.nResidue - 1];
    }

    blobs.clear();  // release blob memory after consuming
}

int Foldcomp::read(std::istream & file, BatchedDecompressionData& batchedData) {
    if (batchedData.is_full()) {
        return -3;
    }

    // Read magic number
    char mNum[MAGICNUMBER_LENGTH];
    file.read(mNum, MAGICNUMBER_LENGTH);
    if (!file) return -1;
    for (int i = 0; i < MAGICNUMBER_LENGTH; i++) {
        if (mNum[i] != MAGICNUMBER[i]) return -1;
    }

    // Read header
    CompressedFileHeader headerLocal;
    file.read(reinterpret_cast<char*>(&headerLocal), sizeof(headerLocal));
    if (!file) return -2;
    if (headerLocal.nResidue <= 0) return -2;

    // nAnchorGroups = inner anchor groups (nAnchor-2) + 1 lastAtomGroup, minimum 1.
    // File anchor coord bytes = nAnchorGroups * 9 floats, all sequential.
    const int nAnchorGroups = (headerLocal.nAnchor <= 1) ? 1 : (headerLocal.nAnchor - 1);

    // Pre-compute payload layout — single allocation below.
    const size_t sz_bb  = (size_t)headerLocal.nResidue          * sizeof(BackboneChain);
    const size_t sz_sc  = (size_t)headerLocal.nSideChainTorsion * sizeof(uint8_t);
    const size_t sz_tf  = (size_t)headerLocal.nResidue          * sizeof(uint8_t);
    const size_t sz_ai  = (size_t)headerLocal.nAnchor           * sizeof(int32_t);
    const size_t sz_ac  = (size_t)nAnchorGroups                 * 9 * sizeof(float);
    const size_t sz_pa  = 9 * sizeof(float);
    const size_t sz_oxt = 3 * sizeof(float);
    const size_t sz_res = (size_t)headerLocal.nResidue          * sizeof(char);
    const size_t sz_rt  = (size_t)headerLocal.nResidue          * sizeof(int8_t);
    const size_t sz_ao  = (size_t)headerLocal.nResidue          * sizeof(int32_t);

    PerStructureBlob blob;
    blob.off_backbone      = 0;
    blob.off_sc_angles     = blob.off_backbone      + sz_bb;
    blob.off_temp_factors  = blob.off_sc_angles     + sz_sc;
    blob.off_anchor_idx    = blob.off_temp_factors  + sz_tf;
    blob.off_anchor_coords = blob.off_anchor_idx    + sz_ai;
    blob.off_prev_atoms    = blob.off_anchor_coords + sz_ac;
    blob.off_oxt_coords    = blob.off_prev_atoms    + sz_pa;
    blob.off_residues      = blob.off_oxt_coords    + sz_oxt;
    blob.off_res_types     = blob.off_residues      + sz_res;
    blob.off_atom_offsets  = blob.off_res_types     + sz_rt;
    const size_t total_payload = blob.off_atom_offsets + sz_ao;

    // Single allocation for all per-file data.
    blob.payload.resize(total_payload);

    // Read anchor indices directly into payload (file order: before title)
    file.read(reinterpret_cast<char*>(blob.payload.data() + blob.off_anchor_idx), sz_ai);
    if (!file) return -2;

    // Read title
    std::string title(headerLocal.lenTitle, '\0');
    if (headerLocal.lenTitle > 0) {
        file.read(&title[0], headerLocal.lenTitle);
        if (!file) return -2;
    }

    // Read prevAtoms directly into payload
    file.read(reinterpret_cast<char*>(blob.payload.data() + blob.off_prev_atoms), sz_pa);
    if (!file) return -2;

    // Read all anchor coordinates in one shot (inner groups + lastAtomGroup, all sequential)
    file.read(reinterpret_cast<char*>(blob.payload.data() + blob.off_anchor_coords), sz_ac);
    if (!file) return -2;

    // Read hasOXT + OXT coordinates directly into payload
    file.read(reinterpret_cast<char*>(&blob.hasOXT), sizeof(char));
    if (!file) return -2;
    file.read(reinterpret_cast<char*>(blob.payload.data() + blob.off_oxt_coords), sz_oxt);
    if (!file) return -2;

    // Read backbone bytes directly into payload
    file.read(reinterpret_cast<char*>(blob.payload.data() + blob.off_backbone), sz_bb);
    if (!file) return -2;

    // Read sidechain angles directly into payload
    if (sz_sc > 0) {
        file.read(reinterpret_cast<char*>(blob.payload.data() + blob.off_sc_angles), sz_sc);
        if (!file) return -2;
    }

    // Read TF calibration scalars + TF bytes directly into payload
    file.read(reinterpret_cast<char*>(&blob.tempFactorDisc_min),  sizeof(float));
    file.read(reinterpret_cast<char*>(&blob.tempFactorDisc_cont), sizeof(float));
    if (!file) return -2;
    if (sz_tf > 0) {
        file.read(reinterpret_cast<char*>(blob.payload.data() + blob.off_temp_factors), sz_tf);
        if (!file) return -2;
    }

    // Fill per-blob metadata
    blob.header            = headerLocal;
    blob.nResidue          = headerLocal.nResidue;
    blob.nAtom             = headerLocal.nAtom;
    blob.nSideChainTorsion = headerLocal.nSideChainTorsion;
    blob.nAllAnchor        = headerLocal.nAnchor;
    blob.nAnchorGroups     = nAnchorGroups;
    blob.firstResidue      = headerLocal.firstResidue;
    blob.lastResidue       = headerLocal.lastResidue;
    blob.useAltAtomOrder   = this->useAltAtomOrder;
    blob.title             = std::move(title);

    // Derive residues, residue types, and per-residue atom-count prefix sums
    // from backbone bytes already in payload — one pass, no extra allocation.
    {
        const uint8_t* bb_bytes = blob.payload.data() + blob.off_backbone;
        char*    res_ptr = reinterpret_cast<char*>(blob.payload.data()   + blob.off_residues);
        int8_t*  rt_ptr  = reinterpret_cast<int8_t*>(blob.payload.data() + blob.off_res_types);
        int32_t* ao_ptr  = reinterpret_cast<int32_t*>(blob.payload.data()+ blob.off_atom_offsets);
        size_t running_atom = 0;
        for (int i = 0; i < headerLocal.nResidue; i++) {
            unsigned int ri = (unsigned int)(bb_bytes[i * 8] >> 3) & 0x1Fu;
            res_ptr[i] = convertIntToOneLetterCode(ri);
            rt_ptr[i]  = static_cast<int8_t>(ri < (unsigned int)GPU_NUM_AMINO_ACIDS
                                              ? ri : (unsigned int)AA_UNK);
            running_atom += getAtomCountFromInt(ri);
            ao_ptr[i] = static_cast<int32_t>(running_atom);
        }
    }

    batchedData.append_blob(std::move(blob));
    return 0;
}

void FoldcompGPUState::initOwnedStreams() {
    createGPUStreamEventPipeline(pipe);
    pipe_owns_streams = true;
}

void FoldcompGPUState::destroyOwnedStreams() {
    destroyGPUStreamEventPipeline(pipe);
    if (e_tf_ready) { cudaEventDestroy(e_tf_ready); e_tf_ready = nullptr; }
    pipe_owns_streams = false;
}

int Foldcomp::decompress_batch_async(AtomCoordinate* atomCoordinates, size_t n_atoms,
                                      cudaEvent_t cpu_sync_event, bool device_emit,
                                      bool emit_pdb, char* pdb_lines_out,
                                      std::vector<int>* out_atom_counts) {
    int success = 0;
    BatchedDecompressionData& batch = gpu.batchData;

    if (batch.current_batch_count == 0) {
        return success;
    }
    gpu.gpu_disc_ctx.setUseFused(gpu.useFusedDisc);

    // Continuized angle pointers — all kept in GPU/pinned memory after continuize_all.
    PinnedFloatVec angles_cont_flat;
    float* tempFactors_cont_flat = nullptr;  // pinned ptr (D→H async on stream A)
    float* d_sc_angles = nullptr;            // device ptr (no D→H at all)
    float* d_tf        = nullptr;            // device ptr (no H→D re-upload to fill kernel)
    size_t n_tf = batch.tempFactorsDiscretized_flat.size();
    float* d_bb_angles = nullptr;
    {
        nvtx3::scoped_range r{"Continuize all angles batch GPU"};
        int total_residues = static_cast<int>(batch.compressedBackBone_flat.size());

        // Build per-struct TF and backbone offsets (1 value per residue each)
        std::vector<int> tf_offsets(batch.current_batch_count + 1, 0);
        for (int i = 0; i < batch.current_batch_count; i++)
            tf_offsets[i + 1] = tf_offsets[i] + batch.nResidue_batch[i];
        // backbone_offsets already contains cumulative residue counts (same as tf_offsets)

        // Build per-struct backbone calibration arrays (6 values per struct)
        std::vector<float> bb_mins_all, bb_cont_fs_all;
        bb_mins_all.reserve(batch.current_batch_count * 6);
        bb_cont_fs_all.reserve(batch.current_batch_count * 6);
        for (int si = 0; si < batch.current_batch_count; si++) {
            const auto& h = batch.headers_batch[si];
            bb_mins_all.insert(bb_mins_all.end(), h.mins, h.mins + 6);
            bb_cont_fs_all.insert(bb_cont_fs_all.end(), h.cont_fs, h.cont_fs + 6);
        }

        FixedAngleDiscretizer scDisc(pow(2, NUM_BITS_TEMP) - 1);

        gpu.gpu_disc_ctx.continuize_all(
            batch.compressedBackBone_flat.data(),  // raw FCZ file bytes (8 bytes/residue)
            total_residues,
            bb_mins_all, bb_cont_fs_all,
            batch.backbone_offsets,      // cumulative residue counts [0, n0, n0+n1, ...]
            angles_cont_flat,            // unused (out_d_bb set)
            batch.tempFactorsDiscretized_flat,
            tf_offsets,
            batch.tempFactorDisc_min_batch, batch.tempFactorDisc_cont_batch,
            batch.sideChainAnglesDiscretized_flat,
            scDisc.min, scDisc.cont_f,
            gpu.streamA(),
            &d_bb_angles, &tempFactors_cont_flat, &d_sc_angles, &d_tf);

        // Record a lightweight event on stream A so we can later sync just up to
        // the TF D→H (which was enqueued by continuize_all), without blocking on
        // the entire stream.  We sync this event before reading tempFactors_cont_flat.
        if (!gpu.e_tf_ready)
            CUDA_CHECK(cudaEventCreateWithFlags(&gpu.e_tf_ready, cudaEventDisableTiming));
        CUDA_CHECK(cudaEventRecord(gpu.e_tf_ready, gpu.streamA()));
    }


    {
        nvtx3::scoped_range r{"GPU Backbone Reconstruction"};

        // Pre-compute total residues for single-pass reservation
        int total_residues = 0;
        for (int i = 0; i < batch.current_batch_count; i++)
            total_residues += batch.nResidue_batch[i];

        // AA type lookup arrays are pre-built by read()/append_batch() on the reader thread.
        // Use batch fields directly — no per-residue CPU work needed here.
        const std::vector<int>&   n_residues_per_protein  = batch.nResidue_batch;
        const std::vector<int8_t>& all_residue_types      = batch.residue_types_flat;
        const std::vector<int>&   all_residue_type_offsets = batch.residue_type_offsets;
        const std::vector<int>&   all_atom_offsets        = batch.atom_offsets_flat;

        // Call GPU backbone reconstruction
        std::vector<int> coord_offsets;

        reconstructBackboneGPU_segmented(
                gpu.gpu_nerf_ctx,
                batch.current_batch_count,
                n_residues_per_protein,
                batch.prevAtoms_flat,
                batch.backbone_offsets,
                all_residue_types,
                all_residue_type_offsets,
                batch.anchorIndices_flat,
                batch.anchor_offsets,
                batch.anchorCoordinates_flat,
                batch.anchor_coord_offsets,
                batch.nAllAnchor_batch,
                gpu.gpu_bb_output_coords,
                coord_offsets,
                d_bb_angles,
                gpu.streamA()
            );

        // Upload backbone metadata to device; sidechain kernels look up backbone inline.
        BBSidechainPtrs bb_ptrs = buildBackboneForSidechains(
            gpu.gpu_nerf_ctx,
            batch.current_batch_count, coord_offsets, n_residues_per_protein,
            batch.prevAtoms_flat, total_residues,
            gpu.streamA()
        );

        // Build per-struct use_alt_order (one bool per struct, not one int per residue)
        std::vector<int8_t> use_alt_per_struct(batch.current_batch_count);
        for (int si = 0; si < batch.current_batch_count; si++)
            use_alt_per_struct[si] = static_cast<int8_t>(batch.useAltAtomOrder_batch[si] ? 1 : 0);

        // Build sidechain metadata on GPU from already-on-device residue types.
        // Computes sc_angle_offsets, output_offsets (prefix sums), and bucket indices.
        SCMetaPtrs sc_meta;
        {
            nvtx3::scoped_range r2{"buildSidechainMeta"};
            sc_meta = buildSidechainMeta(gpu.gpu_sc_ctx, bb_ptrs.d_rt, total_residues, gpu.streamA());
        }

        // GPU Sidechain Reconstruction: all metadata now on device.
        {
            nvtx3::scoped_range r2{"GPU Sidechain Reconstruction Kernel"};
            reconstructSidechainsGPU(
                gpu.gpu_sc_ctx,
                total_residues, batch.current_batch_count,
                bb_ptrs.d_rt,
                bb_ptrs.d_bb_out, bb_ptrs.d_pa, bb_ptrs.d_rfo, bb_ptrs.d_bbo,
                d_sc_angles,
                use_alt_per_struct,
                sc_meta,
                gpu.streamA()
            );
        }

        // Step 3: Fill AtomCoordinate array on GPU
        // all_atom_offsets was built in "bb recon" loop — same data as gpu_output_offsets.
        int total_atoms = all_atom_offsets.back();
        int total_output_atoms = static_cast<int>(n_atoms);

        // Per-struct metadata for the fill kernel.
        std::vector<int>  struct_res_offsets(batch.current_batch_count + 1);
        std::vector<int>  base_anum_per_struct(batch.current_batch_count);
        std::vector<int>  base_ridx_per_struct(batch.current_batch_count);
        std::vector<int>  base_delta_per_struct(batch.current_batch_count);
        // CHAIN_ID_LENGTH bytes/struct (chain_per_struct[si*CHAIN_ID_LENGTH + 0..len-1]),
        // mirroring ChainId (FixedStr<CHAIN_ID_LENGTH>); short chains are zero-padded via
        // the initial '\0' fill below.
        std::vector<char> chain_per_struct(
            static_cast<size_t>(batch.current_batch_count) * CHAIN_ID_LENGTH, '\0');
        struct_res_offsets[0] = 0;
        int global_r  = 0;
        int ac_offset = 0;
        if (out_atom_counts) out_atom_counts->resize(batch.current_batch_count);
        for (int si = 0; si < batch.current_batch_count; si++) {
            int n_res         = batch.nResidue_batch[si];
            int gpu_off_start = all_atom_offsets[global_r];
            base_anum_per_struct[si]  = batch.headers_batch[si].idxAtom;
            base_ridx_per_struct[si]  = batch.headers_batch[si].idxResidue;
            base_delta_per_struct[si] = ac_offset - gpu_off_start;
            {
                const std::string chainName = getChainName(batch.headers_batch[si]);
                for (size_t k = 0; k < CHAIN_ID_LENGTH && k < chainName.size(); k++)
                    chain_per_struct[static_cast<size_t>(si) * CHAIN_ID_LENGTH + k] = chainName[k];
            }
            struct_res_offsets[si + 1] = global_r + n_res;
            int struct_gpu_atoms = all_atom_offsets[global_r + n_res] - gpu_off_start;
            int struct_total_atoms = struct_gpu_atoms + (batch.hasOXT_batch[si] ? 1 : 0);
            // The actual number of atoms this call writes for struct si — may differ
            // from batch.nAtom_batch[si] (the encoded FCZ header's nAtom) if the
            // GPU's own per-residue atom-count table disagrees with the file (e.g. a
            // stale/foreign-encoder file). Callers must slice output using this, not
            // the header value, or every following structure in the batch misaligns.
            if (out_atom_counts) (*out_atom_counts)[si] = struct_total_atoms;
            ac_offset += struct_total_atoms;
            global_r  += n_res;
        }
        // ac_offset is the number of atoms the fill/format kernels below will
        // actually write. It can legitimately fall short of total_output_atoms
        // (e.g. a residue whose GPU-side reconstructed atom count doesn't match
        // the encoded nAtom), which just leaves unused tail capacity — safe.
        // It must never exceed total_output_atoms, or the fill/format kernels
        // would write past the caller-supplied buffer.
        if (ac_offset > total_output_atoms) {
            throw std::runtime_error(
                "decompress_batch_async: computed output atom count (" +
                std::to_string(ac_offset) + ") exceeds caller-supplied n_atoms (" +
                std::to_string(total_output_atoms) + ")");
        }
        if (emit_pdb && pdb_lines_out == nullptr) {
            throw std::runtime_error("decompress_batch_async: pdb_lines_out is null");
        }

        // Sync the TF event before reading tempFactors_cont_flat.
        // By now backbone + sidechain reconstruction (~1ms of GPU work) has been enqueued,
        // so the tiny TF D→H from continuize_all has almost certainly completed.
        CUDA_CHECK(cudaEventSynchronize(gpu.e_tf_ready));

        // Save per-struct OXT state for finish_pending_oxt().
        // ac_offset_base[si] is the cumulative ac_offset at the START of struct si
        // (including GPU atoms + OXT slots from all prior structs).
        {
            auto& p = gpu.pending_oxt;
            p.n_structs = batch.current_batch_count;
            int gr = 0, ao = 0;
            size_t tfo = 0;
            for (int si = 0; si < batch.current_batch_count; si++) {
                int n_res = batch.nResidue_batch[si];
                int sga   = all_atom_offsets[gr + n_res] - all_atom_offsets[gr];
                bool hoxt = static_cast<bool>(batch.hasOXT_batch[si]);
                p.struct_gpu_atoms[si] = sga;
                p.has_oxt[si]          = hoxt;
                p.oxt_xyz[si*3+0]      = batch.OXT_coords_flat[si*3+0];
                p.oxt_xyz[si*3+1]      = batch.OXT_coords_flat[si*3+1];
                p.oxt_xyz[si*3+2]      = batch.OXT_coords_flat[si*3+2];
                {
                    const std::string chainName = getChainName(batch.headers_batch[si]);
                    size_t k = 0;
                    for (; k < CHAIN_ID_LENGTH && k < chainName.size(); k++) p.chain[si][k] = chainName[k];
                    for (; k < CHAIN_ID_LENGTH + 1; k++) p.chain[si][k] = '\0';
                }
                p.last_residue[si]     = batch.headers_batch[si].lastResidue;
                p.idx_atom[si]         = batch.headers_batch[si].idxAtom;
                p.idx_residue[si]      = batch.headers_batch[si].idxResidue;
                p.n_residue[si]        = n_res;
                tfo += n_res;
                float lastTF = (tfo > 0 && tfo - 1 < n_tf)
                                ? tempFactors_cont_flat[tfo - 1] : 0.0f;
                p.last_tf[si]          = lastTF;
                p.ac_offset_base[si]   = ao;
                ao += sga + (hoxt ? 1 : 0);  // advance past GPU atoms + optional OXT slot
                gr += n_res;
            }
        }

        if (emit_pdb) {
            // GPU PDB writer: build d_ac_output on the device, then format the
            // fixed-width ATOM lines on the GPU and D→H only the text (no
            // AtomCoordinate copy). OXT slots are patched host-side at drain.
            nvtx3::scoped_range r2{"GPU PDB writer"};
            AtomCoordinate* d_ac =
                gpu.gpu_sc_ctx.d_ac_output.ensure_n<AtomCoordinate>(total_output_atoms, gpu.streamA());
            // Zero so the OXT slots (not written by the fill kernel) format to
            // blank lines and no garbage FixedStr is scanned out of bounds.
            CUDA_CHECK(cudaMemsetAsync(d_ac, 0,
                (size_t)total_output_atoms * sizeof(AtomCoordinate), gpu.streamA()));
            fillAtomCoordinatesGPU(
                gpu.gpu_sc_ctx,
                total_residues, batch.current_batch_count,
                total_atoms, total_output_atoms,
                static_cast<int>(n_tf), d_tf,
                struct_res_offsets, base_anum_per_struct,
                base_ridx_per_struct, base_delta_per_struct, chain_per_struct,
                nullptr,
                gpu.streamA(), gpu.streamB(),
                gpu.fillDoneEvent(), gpu.slotReadyEvent(),
                cpu_sync_event, /*skip_download=*/true
            );
            char* d_lines = gpu.gpu_pdb_lines_dev.ensure_n<char>(
                (size_t)total_output_atoms * PDB_ATOM_LINE_LEN, gpu.streamA());
            formatAtomLinesDeviceGPU(d_ac, total_output_atoms, d_lines,
                                     gpu.streamA());
            // D→H the formatted lines directly into the caller's shared line
            // buffer at its offset (no per-slot staging + drain capture copy).
            CUDA_CHECK(cudaMemcpyAsync(pdb_lines_out, d_lines,
                (size_t)total_output_atoms * PDB_ATOM_LINE_LEN,
                cudaMemcpyDeviceToHost, gpu.streamA()));
            // Caller syncs cpu_sync_event, then finish_pending_oxt_lines(pdb_lines_out).
            CUDA_CHECK(cudaEventRecord(cpu_sync_event, gpu.streamA()));
        } else if (!device_emit) {
            nvtx3::scoped_range r2{"GPU Output to AtomCoordinate"};
            fillAtomCoordinatesGPU(
                gpu.gpu_sc_ctx,
                total_residues, batch.current_batch_count,
                total_atoms, total_output_atoms,
                static_cast<int>(n_tf), d_tf,
                struct_res_offsets, base_anum_per_struct,
                base_ridx_per_struct, base_delta_per_struct, chain_per_struct,
                atomCoordinates,
                gpu.streamA(), gpu.streamB(),
                gpu.fillDoneEvent(), gpu.slotReadyEvent(),
                cpu_sync_event
            );
            // D→H is now running on stream B.  Caller must wait on cpu_sync_event
            // (cudaEventSynchronize) then call finish_pending_oxt() before reading
            // atomCoordinates.
        } else {
            // Device-emit path: pack xyz into gpu_packed_coords in final output
            // order (incl OXT), keeping coordinates resident on the device.
            // reconstructSidechainsGPU already wrote gpu_sc_ctx.output_coords
            // (float[total_atoms*3], GPU-atom order) on streamA; each structure's
            // atoms are contiguous there, so a per-structure D→D copy suffices.
            nvtx3::scoped_range r2{"GPU Pack Coords (device emit)"};
            float* d_packed =
                gpu.gpu_packed_coords.ensure_n<float>((size_t)total_output_atoms * 3, gpu.streamA());
            const float* d_src =
                static_cast<const float*>(gpu.gpu_sc_ctx.output_coords.ptr);
            for (int si = 0; si < batch.current_batch_count; si++) {
                const int gpu_off = all_atom_offsets[struct_res_offsets[si]];
                const int dst_off = gpu.pending_oxt.ac_offset_base[si];
                const int n_gpu   = gpu.pending_oxt.struct_gpu_atoms[si];
                if (n_gpu > 0) {
                    CUDA_CHECK(cudaMemcpyAsync(
                        d_packed + (size_t)dst_off * 3,
                        d_src + (size_t)gpu_off * 3,
                        (size_t)n_gpu * 3 * sizeof(float),
                        cudaMemcpyDeviceToDevice, gpu.streamA()));
                }
                if (gpu.pending_oxt.has_oxt[si]) {
                    CUDA_CHECK(cudaMemcpyAsync(
                        d_packed + (size_t)(dst_off + n_gpu) * 3,
                        &gpu.pending_oxt.oxt_xyz[si * 3],
                        3 * sizeof(float),
                        cudaMemcpyHostToDevice, gpu.streamA()));
                }
            }
            // Caller syncs cpu_sync_event before reading gpu_packed_coords.
            CUDA_CHECK(cudaEventRecord(cpu_sync_event, gpu.streamA()));
        }
    }

    return success;
}

void Foldcomp::finish_pending_oxt(AtomCoordinate* atomCoordinates) {
    // cpu_sync_event must have been synchronized by the caller before invoking this.
    const auto& p = gpu.pending_oxt;
    nvtx3::scoped_range r2{"Finishing pending OXT"};
    for (int si = 0; si < p.n_structs; si++) {
        if (!p.has_oxt[si]) continue;
        // OXT slot is at ac_offset_base[si] + sga (immediately after this struct's GPU atoms).
        int ac_off = p.ac_offset_base[si] + p.struct_gpu_atoms[si];
        std::string lastResStr = getThreeLetterCode(p.last_residue[si]);
        int oxt_atom_num = p.idx_atom[si] + p.struct_gpu_atoms[si];
        int oxt_residue_num = p.idx_residue[si] + p.n_residue[si] - 1;
        AtomCoordinate oxtAtom("OXT", lastResStr, p.chain[si], oxt_atom_num, oxt_residue_num,
                                p.oxt_xyz[si*3+0], p.oxt_xyz[si*3+1], p.oxt_xyz[si*3+2]);
        oxtAtom.tempFactor = p.last_tf[si];
        atomCoordinates[ac_off] = oxtAtom;
    }
}

void Foldcomp::finish_pending_oxt_lines(char* lines) {
    // cpu_sync_event must have been synchronized by the caller before invoking this.
    const auto& p = gpu.pending_oxt;
    for (int si = 0; si < p.n_structs; si++) {
        if (!p.has_oxt[si]) continue;
        int ac_off = p.ac_offset_base[si] + p.struct_gpu_atoms[si];
        std::string lastResStr = getThreeLetterCode(p.last_residue[si]);
        int oxt_atom_num = p.idx_atom[si] + p.struct_gpu_atoms[si];
        int oxt_residue_num = p.idx_residue[si] + p.n_residue[si] - 1;
        AtomCoordinate oxtAtom("OXT", lastResStr, p.chain[si], oxt_atom_num, oxt_residue_num,
                                p.oxt_xyz[si*3+0], p.oxt_xyz[si*3+1], p.oxt_xyz[si*3+2]);
        oxtAtom.tempFactor = p.last_tf[si];
        formatAtomLineHost(oxtAtom, lines + (size_t)ac_off * PDB_ATOM_LINE_LEN);
    }
}

int Foldcomp::decompress_batch(AtomCoordinate* atomCoordinates, size_t n_atoms) {
    cudaEvent_t cpu_evt = acquireCPUSyncEvent();
    int r = decompress_batch_async(atomCoordinates, n_atoms, cpu_evt);
    CUDA_CHECK(cudaEventSynchronize(cpu_evt));
    releaseCPUSyncEvent(cpu_evt);
    finish_pending_oxt(atomCoordinates);
    return r;
}
