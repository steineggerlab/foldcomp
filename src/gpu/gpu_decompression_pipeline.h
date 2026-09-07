// SPDX-License-Identifier: MIT

#pragma once

#ifdef FOLDCOMP_WITH_CUDA

#include "structure_codec.h"

#include <functional>
#include <string>
#include <vector>

class Processor;

struct GPUDecompressionConfig {
    int gpu_slots = 3;
    int reader_threads = 4;
    int max_batch_size = 100;
    int max_residues = 3000;
    bool use_alt_order = false;
    bool use_fused_discretizer = true;
    // When true, the pipeline formats PDB ATOM lines on the GPU and delivers them
    // (plus per-segment metadata) via GPUDecompressionResult::pdb_lines /
    // line_segments instead of host DecompressedSegments. Uses the host `output`
    // callback (segments stay empty).
    bool emit_device_pdb = false;
    // Max number of distinct (gpu_slots, reader_threads, max_batch_size, max_residues,
    // use_fused_discretizer, mode, device) pipeline instances kept resident at once.
    // Concurrent runGPUDecompressionPipeline calls with distinct configs (or the same
    // config, if more than pool_size calls are in flight) genuinely overlap on the GPU
    // up to this many at a time; beyond that, callers block until an instance frees up.
    // Each resident instance duplicates the full GPU-slot + batch-pool preallocation, so
    // raising this trades VRAM/pinned-memory footprint for concurrency.
    int pool_size = 2;

    // Validates gpu_slots/reader_threads/max_batch_size/max_residues/pool_size are
    // positive and max_batch_size fits FoldcompGPUState::MAX_STRUCTS_PER_BATCH. Every
    // caller that builds a GPUDecompressionConfig should call this instead of
    // duplicating these checks, so validation stays consistent in one place.
    // Defined in gpu_decompression_pipeline.cpp (needs FoldcompGPUState::MAX_STRUCTS_PER_BATCH).
    bool validate(std::string& error) const;
};

// A fragment's location within a pipeline instance's PDB line capture buffer
// (offset and length in atoms; each atom is one PDB_ATOM_LINE_LEN-byte line).
struct GPUCoordFragment {
    size_t offset = 0;
    int    n_atoms = 0;
};

struct GPUDecompressionResult {
    std::string name;
    std::string title;
    std::vector<DecompressedSegment> segments;
    // Populated only in emit_device_pdb mode: per-segment metadata plus, for each
    // fragment (parallel to line_segments), its offset into the pipeline instance's
    // line capture buffer (GPUPipelineLease::pdbCaptureBase). The stitch reads
    // lines from there.
    std::vector<PdbLineSegment> line_segments;
    std::vector<GPUCoordFragment> pdb_frags;
};

using GPUDecompressionOutputFn =
    std::function<bool(GPUDecompressionResult&& result)>;

// Opaque handle to a pooled GPU pipeline instance (see the pool implementation in
// gpu_decompression_pipeline.cpp), returned by runGPUDecompressionPipeline when a
// caller passes a non-null out_lease. While the lease is alive, the instance stays
// checked out of the pool — no other thread can reuse or evict it — so its capture
// buffer (populated by an emit_device_pdb run) is safe to read with no additional
// locking. Destroying the lease (e.g. end of scope) returns the instance to the
// pool. Move-only.
class GPUPipelineLease {
public:
    GPUPipelineLease() = default;
    ~GPUPipelineLease();
    GPUPipelineLease(const GPUPipelineLease&) = delete;
    GPUPipelineLease& operator=(const GPUPipelineLease&) = delete;
    GPUPipelineLease(GPUPipelineLease&& other) noexcept : pipe_(other.pipe_) {
        other.pipe_ = nullptr;
    }
    GPUPipelineLease& operator=(GPUPipelineLease&& other) noexcept;

    // Base pointer of this instance's PDB line capture buffer (PDB_ATOM_LINE_LEN
    // bytes/atom). Only valid while this lease is alive.
    const char* pdbCaptureBase() const;

    // Opaque; set internally by runGPUDecompressionPipelineImpl (points at a
    // PooledPipeline, whose full definition lives only in gpu_decompression_pipeline.cpp).
    void* pipe_ = nullptr;
};

bool cudaRuntimeAvailable(std::string* reason = nullptr);

// Free every pooled GPU pipeline instance (GPU slots + batch pool + capture
// buffers) that runGPUDecompressionPipeline reuses across calls. Waits for any
// in-flight runs to finish first. Deterministically releases the retained
// (pinned + device) buffers. Safe to call between runs; the next run rebuilds
// the pool on demand.
void releaseGPUPipelineCache();

// out_lease, if non-null, receives ownership of the pipeline instance that served
// this run instead of returning it to the pool — use this when the caller needs to
// read GPUDecompressionResult::pdb_frags out of the instance's capture buffer after
// the run returns (see GPUPipelineLease above).
bool runGPUDecompressionPipeline(
    Processor& processor,
    const GPUDecompressionConfig& config,
    const GPUDecompressionOutputFn& output,
    std::string& error,
    GPUPipelineLease* out_lease = nullptr);

#endif
