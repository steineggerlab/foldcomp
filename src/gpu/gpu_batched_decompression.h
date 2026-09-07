// SPDX-License-Identifier: MIT
/**
 * File: gpu_batched_decompression.h
 * Project: foldcomp
 * Description:
 *     GPU batch data structures used by the batched decompression pipeline.
 *     Extracted from foldcomp.h. This header is only ever included from
 *     foldcomp.h, at the point where CompressedFileHeader, BackboneChain,
 *     AminoAcidIndex and CudaPinnedAllocator (from gpu_buffer.h) are already
 *     visible in the translation unit — it does not re-include foldcomp.h
 *     itself to avoid a circular include.
 *
 *     Included unconditionally (not just under FOLDCOMP_WITH_CUDA) so that
 *     Foldcomp's batch-decompression method declarations are identical across
 *     CUDA and CPU-only builds: BatchedDecompressionData is forward-declared
 *     here for both, and CPU-only builds get a cudaEvent_t stand-in plus a
 *     throwGpuBatchUnsupported() helper used by the CPU-only method stubs
 *     declared in foldcomp.h.
 */
#pragma once

// Forward declaration only; the full definition (below) is CUDA-only. Lets
// Foldcomp::read(istream&, BatchedDecompressionData&) keep the same signature
// in CPU-only builds, where it's never actually called.
struct BatchedDecompressionData;

#ifndef FOLDCOMP_WITH_CUDA
#include <stdexcept>
#include <string>

// Stand-in for CUDA's cudaEvent_t (itself a pointer-to-incomplete-struct
// typedef), so Foldcomp::decompress_batch_async keeps an identical signature
// without pulling in CUDA headers.
using cudaEvent_t = void*;

// Thrown by the CPU-only stub implementations of Foldcomp's GPU batch-
// decompression methods (decompress_batch, etc.) declared in foldcomp.h.
// Those stubs exist purely so Foldcomp's public API is identical between
// CUDA and CPU-only builds (MR !1 comment c0c03bb7); none of them do
// anything without CUDA.
[[noreturn]] inline void throwGpuBatchUnsupported(const char* method) {
    throw std::runtime_error(
        std::string("Foldcomp::") + method + " requires a CUDA-enabled build (FOLDCOMP_WITH_CUDA)");
}
#endif

#ifdef FOLDCOMP_WITH_CUDA
/**
 * All data for a single parsed FCZ structure packed into one contiguous
 * allocation. Produced by Foldcomp::read(istream, BatchedDecompressionData)
 * on reader threads; consumed by BatchedDecompressionData::finalize() on the
 * GPU consumer thread to build the flat GPU-upload arrays in bulk.
 *
 * Payload layout (byte offsets computed at construction time):
 *   [off_backbone]      BackboneChain[nResidue]        (8 bytes each)
 *   [off_sc_angles]     uint8_t[nSideChainTorsion]
 *   [off_temp_factors]  uint8_t[nResidue]
 *   [off_anchor_idx]    int32_t[nAllAnchor]
 *   [off_anchor_coords] float[nAllAnchor * 3]
 *   [off_prev_atoms]    float[9]   (always 36 bytes)
 *   [off_oxt_coords]    float[3]   (always 12 bytes)
 *   [off_residues]      char[nResidue]    (one-letter codes)
 *   [off_res_types]     int8_t[nResidue] (AminoAcidIndex)
 *   [off_atom_offsets]  int32_t[nResidue] (local cumulative atom counts from 0)
 */
struct PerStructureBlob {
    CompressedFileHeader header;
    int nResidue = 0, nAtom = 0, nSideChainTorsion = 0, nAllAnchor = 0;
    int nAnchorGroups = 0;  // = anchorCoordinates.size() = nAllAnchor-1 (groups of 9 floats)
    int firstResidue = 0, lastResidue = 0;
    float tempFactorDisc_min = 0.f, tempFactorDisc_cont = 0.f;
    bool hasOXT = false, useAltAtomOrder = false;
    std::string title;

    std::vector<uint8_t> payload;
    size_t off_backbone = 0, off_sc_angles = 0, off_temp_factors = 0;
    size_t off_anchor_idx = 0, off_anchor_coords = 0;
    size_t off_prev_atoms = 0, off_oxt_coords = 0;
    size_t off_residues = 0, off_res_types = 0, off_atom_offsets = 0;

    // backbone_raw() returns raw FCZ file bytes (8 bytes/residue, cross-byte bit layout).
    // NOT struct-format BackboneChain — passed directly to the GPU per-struct unpack kernel.
    const uint64_t*      backbone_raw()  const { return reinterpret_cast<const uint64_t*>(payload.data() + off_backbone); }
    const uint8_t*       sc_angles()     const { return payload.data() + off_sc_angles; }
    const uint8_t*       temp_factors()  const { return payload.data() + off_temp_factors; }
    const int32_t*       anchor_idx()    const { return reinterpret_cast<const int32_t*>(payload.data() + off_anchor_idx); }
    const float*         anchor_coords() const { return reinterpret_cast<const float*>(payload.data() + off_anchor_coords); }
    const float*         prev_atoms()    const { return reinterpret_cast<const float*>(payload.data() + off_prev_atoms); }
    const float*         oxt_coords()    const { return reinterpret_cast<const float*>(payload.data() + off_oxt_coords); }
    const char*          residues()      const { return reinterpret_cast<const char*>(payload.data() + off_residues); }
    const int8_t*        res_types()     const { return reinterpret_cast<const int8_t*>(payload.data() + off_res_types); }
    const int32_t*       atom_offsets()  const { return reinterpret_cast<const int32_t*>(payload.data() + off_atom_offsets); }
};

// Batch structures for efficient batched decompression
// All data is stored in flat arrays with offsets for each structure in the batch
struct BatchedDecompressionData {
    // Batch configuration
    int max_batch_size = 100;
    int current_batch_count{0};

    BatchedDecompressionData() {
        clear();
    }

    // Offsets into flat arrays for each structure in batch
    // Size: max_batch_size + 1 (last element is total size)
    std::vector<int> backbone_offsets;          // Offset into compressedBackBone_flat (angle offsets = *6)
    std::vector<int> sidechain_offsets;         // Offset into sideChainAnglesDiscretized_flat
    std::vector<int> anchor_offsets;            // Offset into anchorIndices_flat
    std::vector<int> anchor_coord_offsets;      // Offset into anchorCoordinates_flat

    // Per-structure metadata (size: max_batch_size)
    std::vector<int> nResidue_batch;
    std::vector<int> nAtom_batch;
    std::vector<int> nSideChainTorsion_batch;
    std::vector<int> nAllAnchor_batch;
    std::vector<char> firstResidue_batch;
    std::vector<char> lastResidue_batch;
    std::vector<char> hasOXT_batch;
    std::vector<bool> useAltAtomOrder_batch;
    std::vector<std::string> titles_batch;
    std::vector<float> tempFactorDisc_min_batch;
    std::vector<float> tempFactorDisc_cont_batch;

    // Headers (size: max_batch_size)
    std::vector<CompressedFileHeader> headers_batch;

    // Flat compressed backbone data — raw FCZ file bytes (8 bytes/residue, cross-byte bit layout).
    // Stored as uint64_t (8 bytes each) for direct H→D upload; interpreted by
    // unpack_raw_file_continuize_bb_per_struct_kernel in gpu_discretizer.cu.
    // Total size: sum of all nResidue values
    std::vector<uint64_t, CudaPinnedAllocator<uint64_t>> compressedBackBone_flat;

    // Flat sidechain data (uint8_t: values are 4-bit, 0-15)
    // Total size: sum of all nSideChainTorsion values
    std::vector<uint8_t, CudaPinnedAllocator<uint8_t>> sideChainAnglesDiscretized_flat;

    // Flat temperature factor data (uint8_t: values are 8-bit, 0-255)
    // Total size: sum of all nResidue values
    std::vector<uint8_t> tempFactorsDiscretized_flat;

    // Flat anchor indices
    // Total size: sum of all nAllAnchor values
    std::vector<int> anchorIndices_flat;

    // Flat anchor coordinates (3 atoms * 3 coords per anchor point)
    // Total size: sum of all (nAllAnchor * 9) values
    std::vector<float> anchorCoordinates_flat;

    // Flat prev atoms (3 atoms * 3 coords per structure)
    // Size: max_batch_size * 9
    std::vector<float> prevAtoms_flat;

    // Flat OXT coordinates
    // Size: max_batch_size * 3
    std::vector<float> OXT_coords_flat;

    // Residue names (stored as chars for efficiency)
    // Total size: sum of all nResidue values
    std::vector<char> residues_flat;

    // Pre-computed AA lookup data (built during read()/append_batch())
    // Avoids O(total_residues) CPU work in decompress_batch on the GPU dispatch thread.
    // Layout per struct: [BB[0]_aa, ..., BB[nResidue-1]_aa] — nResidue entries
    // compressedBackBone[i].residue encodes the AA type for residue i (including residue 0).
    // firstResidue is NOT pushed separately; doing so would double-count residue 0.
    std::vector<int8_t> residue_types_flat;     // [sum(nResidue)] AminoAcidIndex per residue
    std::vector<int>    residue_type_offsets;   // [n_structs+1] per-struct start in residue_types_flat
    std::vector<int>    atom_offsets_flat;      // [sum(nResidue)+1] prefix sums of ATOM_COUNTS
    // Running prefix sum carried across read() calls
    int    _atom_running{0};

    // Problem 3: pre-finalize blobs (one per parsed structure, moved in O(1)).
    // Flat arrays above are empty until finalize() is called.
    std::vector<PerStructureBlob> blobs;
    int cached_total_atoms{0};  // sum of nAtom across blobs; valid before finalize()

    // Methods for batch management
    void clear() {
        current_batch_count = 0;
        backbone_offsets.clear();
        sidechain_offsets.clear();
        anchor_offsets.clear();
        anchor_coord_offsets.clear();

        nResidue_batch.clear();
        nAtom_batch.clear();
        nSideChainTorsion_batch.clear();
        nAllAnchor_batch.clear();
        firstResidue_batch.clear();
        lastResidue_batch.clear();
        hasOXT_batch.clear();
        useAltAtomOrder_batch.clear();
        titles_batch.clear();
        tempFactorDisc_min_batch.clear();
        tempFactorDisc_cont_batch.clear();

        headers_batch.clear();

        compressedBackBone_flat.clear();
        sideChainAnglesDiscretized_flat.clear();
        tempFactorsDiscretized_flat.clear();
        anchorIndices_flat.clear();
        anchorCoordinates_flat.clear();
        prevAtoms_flat.clear();
        OXT_coords_flat.clear();
        residues_flat.clear();

        residue_types_flat.clear();
        residue_type_offsets.clear();
        atom_offsets_flat.clear();
        _atom_running = 0;
        blobs.clear();
        cached_total_atoms = 0;

        // Reinitialize offset arrays to contain starting zero
        backbone_offsets.push_back(0);
        sidechain_offsets.push_back(0);
        anchor_offsets.push_back(0);
        anchor_coord_offsets.push_back(0);
        residue_type_offsets.push_back(0);
        atom_offsets_flat.push_back(0);
    }

    void reserve() {
        clear();

        // Reserve space for offsets
        backbone_offsets.reserve(max_batch_size + 1);
        sidechain_offsets.reserve(max_batch_size + 1);
        anchor_offsets.reserve(max_batch_size + 1);
        anchor_coord_offsets.reserve(max_batch_size + 1);

        // Reserve space for per-structure metadata
        nResidue_batch.reserve(max_batch_size);
        nAtom_batch.reserve(max_batch_size);
        nSideChainTorsion_batch.reserve(max_batch_size);
        nAllAnchor_batch.reserve(max_batch_size);
        firstResidue_batch.reserve(max_batch_size);
        lastResidue_batch.reserve(max_batch_size);
        hasOXT_batch.reserve(max_batch_size);
        useAltAtomOrder_batch.reserve(max_batch_size);
        titles_batch.reserve(max_batch_size);
        tempFactorDisc_min_batch.reserve(max_batch_size);
        tempFactorDisc_cont_batch.reserve(max_batch_size);
        headers_batch.reserve(max_batch_size);

        // Flat arrays will grow dynamically
        prevAtoms_flat.reserve(max_batch_size * 9);
        OXT_coords_flat.reserve(max_batch_size * 3);
        blobs.reserve(max_batch_size);
    }

    bool is_full() const {
        return current_batch_count >= max_batch_size;
    }

    /**
     * Move a parsed blob into this batch in O(1) — no flat array copies.
     * Flat arrays are populated lazily by finalize() before GPU dispatch.
     */
    void append_blob(PerStructureBlob&& b) {
        cached_total_atoms += b.nAtom;
        blobs.push_back(std::move(b));
        current_batch_count++;
    }

    /**
     * Build all flat GPU-upload arrays from stored blobs using pre-reserved
     * bulk insertions. No-op if blobs is empty (safe to call multiple times).
     * Must be called before decompress_batch().
     */
    void finalize();

    /**
     * Pre-allocate pinned flat arrays to their expected steady-state sizes.
     * Call once per slot at startup to avoid cold-start cudaMallocHost latency
     * inside the first finalize() calls.
     *
     * @param max_total_residues  Upper bound on total residues per batch
     * @param max_total_sc        Upper bound on total sidechain torsion angles per batch
     */
    void preallocate_pinned(size_t max_total_residues, size_t max_total_sc) {
        compressedBackBone_flat.reserve(max_total_residues);
        sideChainAnglesDiscretized_flat.reserve(max_total_sc);
    }
};

/**
 * All per-instance GPU-only state for a Foldcomp object: batch data, grow-only
 * device/pinned buffers, and the per-instance stream/event pipeline. Grouped
 * here (rather than as flat members on Foldcomp) so the CPU-only build sees
 * no trace of it beyond the single `#ifdef FOLDCOMP_WITH_CUDA` `gpu` member
 * on Foldcomp.
 */
struct FoldcompGPUState {
    // Batched decompression data
    BatchedDecompressionData batchData;

    bool isBatched{true};
    bool useFusedDisc = true;         // use fused uint→float+scale kernel for continuization

    // Persistent GPU context: all grow-only device/pinned buffers for this Foldcomp instance.
    GPUDiscretizerContext gpu_disc_ctx;
    NerfBuffers           gpu_nerf_ctx;
    SidechainBuffers      gpu_sc_ctx;
    PinnedFloatVec        gpu_bb_output_coords; // backbone output staging (skip_download path)
    GrowOnlyPinnedBuf     gpu_ac_buf;  // AtomCoordinate output staging

    // Device-emit path: packed xyz (float[total_output_atoms*3]) in final output
    // order (incl OXT), produced by decompress_batch_async when device_emit=true.
    // Grow-only, per-slot; the pipeline copies it out at drain before slot reuse.
    GrowOnlyDeviceBuf     gpu_packed_coords;

    // GPU PDB writer path (emit_pdb=true): PDB ATOM lines (81 bytes/atom) formatted
    // on the device from d_ac_output, then D→H'd into the caller-provided shared
    // line buffer. OXT slots are patched host-side by finish_pending_oxt_lines.
    GrowOnlyDeviceBuf     gpu_pdb_lines_dev;

    // Two-stream pipeline state: OXT data saved by decompress_batch_async.
    static constexpr int MAX_STRUCTS_PER_BATCH = 100; // matches default max_batch_size
    struct PendingOXTState {
        int   n_structs            = 0;
        int   struct_gpu_atoms[MAX_STRUCTS_PER_BATCH];
        bool  has_oxt[MAX_STRUCTS_PER_BATCH];
        float oxt_xyz[MAX_STRUCTS_PER_BATCH * 3]; // flat: [si*3 + 0..2]
        char  chain[MAX_STRUCTS_PER_BATCH][CHAIN_ID_LENGTH + 1]; // mirrors ChainId (FixedStr<CHAIN_ID_LENGTH>)
        char  last_residue[MAX_STRUCTS_PER_BATCH];
        int   idx_atom[MAX_STRUCTS_PER_BATCH];
        int   idx_residue[MAX_STRUCTS_PER_BATCH];
        int   n_residue[MAX_STRUCTS_PER_BATCH];
        float last_tf[MAX_STRUCTS_PER_BATCH];
        int   ac_offset_base[MAX_STRUCTS_PER_BATCH]; // cumul ac_offset at start of struct si
    } pending_oxt;

    cudaEvent_t e_tf_ready = nullptr; // recorded on stream A after TF D→H in continuize_all

    // Per-instance streams and events for multi-pipeline concurrent execution.
    // Call initOwnedStreams() to create them (every FoldcompGPUState built by the
    // pipeline pool does, in buildPooledPipeline); the create/destroy logic is
    // shared via createGPUStreamEventPipeline()/destroyGPUStreamEventPipeline()
    // (gpu_stream.h).
    GPUStreamEventPipeline pipe;
    bool pipe_owns_streams = false;

    // Returns the compute stream for this instance. Valid only after initOwnedStreams().
    cudaStream_t streamA()     const { return pipe.stream_a; }
    // Returns the D→H stream for this instance. Valid only after initOwnedStreams().
    cudaStream_t streamB()     const { return pipe.stream_b; }
    // Returns the fill-done event for this instance. Valid only after initOwnedStreams().
    cudaEvent_t  fillDoneEvent()  const { return pipe.fill_done; }
    // Returns the slot-ready guard event for this instance. Valid only after initOwnedStreams().
    cudaEvent_t  slotReadyEvent() const { return pipe.slot_ready; }

    // Create owned streams and events; pre-record slot_ready as "completed".
    void initOwnedStreams();
    // Destroy owned streams and events (called automatically by the destructor if
    // initOwnedStreams was called).
    void destroyOwnedStreams();

    ~FoldcompGPUState() { if (pipe_owns_streams) destroyOwnedStreams(); }
};
#endif
