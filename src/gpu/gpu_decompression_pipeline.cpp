// SPDX-License-Identifier: MIT

#include "gpu_decompression_pipeline.h"

#ifdef FOLDCOMP_WITH_CUDA

#include "foldcomp.h"
#include "gpu_pdb_writer.h"
#include "gpu_sidechain.h"
#include "gpu_stream.h"
#include "input_processor.h"

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <cstring>
#include <deque>
#include <exception>
#include <memory>
#include <mutex>
#include <queue>
#include <sstream>
#include <thread>

namespace {

struct MemBuf : std::streambuf {
    MemBuf(const char* base, size_t size) {
        char* p = const_cast<char*>(base);
        setg(p, p, p + size);
    }
};

struct MemStream : MemBuf, std::istream {
    MemStream(const char* base, size_t size)
        : MemBuf(base, size), std::istream(static_cast<std::streambuf*>(this)) {}
};

struct LogicalAssembly {
    std::string name;
    std::string title;
    std::vector<DecompressedSegment> slots;
    int pending_gpu = 0;
    bool emitted = false;
    // emit_device_pdb only: per-slot offset into the shared line capture buffer for
    // GPU fragments (SIZE_MAX = raw/CPU-fallback slot, appended at emit) + metadata.
    std::vector<size_t>         slot_pdb_off;
    std::vector<PdbLineSegment> slot_pdb_meta;
};

struct FragmentTarget {
    std::shared_ptr<LogicalAssembly> assembly;
    size_t slot = 0;
    int model = 1;
    std::string chain;
    bool override_chain = false;
};

struct ParsedGPUFragment {
    PerStructureBlob blob;
    FragmentTarget target;
};

struct ParsedInput {
    std::shared_ptr<LogicalAssembly> assembly;
    std::vector<ParsedGPUFragment> gpu_fragments;
};

bool parseFczBlob(
    Foldcomp& parser,
    const char* data,
    size_t size,
    bool use_alt_order,
    PerStructureBlob& blob,
    std::string& error) {
    BatchedDecompressionData parsed;
    parser.useAltAtomOrder = use_alt_order;
    MemStream stream(data, size);
    int flag = parser.read(stream, parsed);
    if (flag != 0 || parsed.blobs.size() != 1) {
        error = "invalid FCZ fragment (read error " + std::to_string(flag) + ")";
        return false;
    }
    blob = std::move(parsed.blobs.front());
    return true;
}

bool parseInput(
    Foldcomp& parser,
    const char* name,
    const char* data,
    size_t size,
    const GPUDecompressionConfig& config,
    ParsedInput& parsed,
    std::string& error) {
    parsed = {};
    parsed.assembly = std::make_shared<LogicalAssembly>();
    parsed.assembly->name = name ? name : "";

    if (!hasContainerMagic(data, size)) {
        PerStructureBlob blob;
        if (!parseFczBlob(parser, data, size, config.use_alt_order, blob, error)) {
            return false;
        }
        parsed.assembly->title = blob.title;
        parsed.assembly->slots.resize(1);
        parsed.assembly->pending_gpu = 1;
        FragmentTarget target;
        target.assembly = parsed.assembly;
        parsed.gpu_fragments.push_back({std::move(blob), std::move(target)});
        return true;
    }

    std::vector<ContainerFragment> fragments;
    if (!readContainer(data, size, parsed.assembly->title, fragments)) {
        error = "invalid FCZC container";
        return false;
    }
    if (fragments.empty()) {
        error = "FCZC container has no fragments";
        return false;
    }

    parsed.assembly->slots.resize(fragments.size());
    parsed.gpu_fragments.reserve(fragments.size());
    for (size_t i = 0; i < fragments.size(); ++i) {
        const ContainerFragment& fragment = fragments[i];
        if (fragment.kind == CONTAINER_FRAGMENT_KIND_RAW_ATOMS) {
            std::vector<AtomCoordinate> atoms;
            if (!deserializeAtomCoordinates(
                    fragment.payload.data(), fragment.payload.size(), atoms)) {
                error = "invalid raw-atom fragment at index " + std::to_string(i);
                return false;
            }
            for (AtomCoordinate& atom : atoms) {
                atom.model = fragment.model;
                atom.chain = fragment.chain;
            }
            parsed.assembly->slots[i] = {fragment.model, std::move(atoms)};
            continue;
        }
        if (fragment.kind != CONTAINER_FRAGMENT_KIND_FCZ) {
            error = "unsupported container fragment kind " +
                    std::to_string(fragment.kind) + " at index " + std::to_string(i);
            return false;
        }

        PerStructureBlob blob;
        if (!parseFczBlob(parser, fragment.payload.data(), fragment.payload.size(),
                config.use_alt_order, blob, error)) {
            error += " at container fragment " + std::to_string(i);
            return false;
        }
        FragmentTarget target;
        target.assembly = parsed.assembly;
        target.slot = i;
        target.model = fragment.model;
        target.chain = fragment.chain;
        target.override_chain = true;
        parsed.gpu_fragments.push_back({std::move(blob), std::move(target)});
        parsed.assembly->pending_gpu++;
    }
    return true;
}

} // namespace

bool cudaRuntimeAvailable(std::string* reason) {
    int count = 0;
    cudaError_t status = cudaGetDeviceCount(&count);
    if (status != cudaSuccess || count == 0) {
        if (reason) {
            *reason = status == cudaSuccess
                          ? "no CUDA devices are available"
                          : cudaGetErrorString(status);
        }
        cudaGetLastError();
        return false;
    }
    if (reason) reason->clear();
    return true;
}

// Batch of parsed structures dispatched to one GPU slot.
struct Batch {
    BatchedDecompressionData data;
    std::vector<FragmentTarget> targets;
};

// Bookkeeping-only mutex: guards the pool container (g_pool) and instance
// checkout/return, never held for the duration of a run — nor, since an instance is
// reserved (in_use = true) before the lock is released, for the duration of a build
// or evict+rebuild either (see buildPooledPipelineUnlocked), so warming up pipelines
// on different devices concurrently doesn't serialize on each other's cudaMallocHost
// preallocate. File scope so releaseGPUPipelineCache() can lock it too.
static std::mutex g_pool_mutex;
static std::condition_variable g_pool_cv;
static uint64_t g_pool_tick = 0;

// Pipeline output mode: host = build host DecompressedSegments (+ optional GPU PDB
// lines), pdb = GPU-formatted PDB lines.
enum PipeMode { MODE_HOST = 0, MODE_PDB = 1 };

// One pooled GPU pipeline instance: GPU slots + batch pool (rebuilding these every
// call is dominated by cudaMallocHost pinned preallocate, so pooling removes that
// per-call setup floor) plus this instance's own output capture buffers, so
// concurrently-running instances never share mutable state. Fields other than
// in_use/last_used_tick/device_id/gpu_slots are only touched by the thread that
// currently holds this instance checked out (see checkoutPipeline/returnPipeline
// below); in_use and last_used_tick are only touched under g_pool_mutex. device_id
// and gpu_slots are set under g_pool_mutex too — before the lock is released for a
// build (see buildPooledPipelineUnlocked) — because growCpuSyncPoolForResidentInstances
// reads them for every pool entry regardless of in_use.
struct PooledPipeline {
    int gpu_slots = 0, reader_threads = 0, max_batch_size = 0, max_residues = 0;
    bool use_fused = true;
    int mode = MODE_HOST;
    int device_id = -1;
    std::vector<std::unique_ptr<Foldcomp>> gpu_fc;
    std::vector<std::unique_ptr<Batch>> batch_pool;

    // Per-instance PINNED host buffer that the GPU PDB writer D→H's all formatted
    // ATOM lines into (grow-with-preserve). The drain records per-fragment offsets
    // instead of copying, and the parallel stitch reads straight from here.
    char*  pdb_capture = nullptr;
    size_t pdb_capacity = 0; // bytes
    size_t pdb_used = 0;     // atoms (lines) written so far this run

    // Per-slot "last append landed" event for the capture buffer above. Recorded on
    // a slot's own compute stream right after decompress_batch_async, whose internal
    // D→H write lands on the same stream. Growth waits on just these (<= MAX_IN_FLIGHT
    // of them) instead of cudaDeviceSynchronize(), so a capture-buffer resize only
    // blocks on the slots that actually write into it, not unrelated device work (see
    // MR !1 bdcb5db2: stream-ordered growth vs. device-wide sync).
    cudaEvent_t pdb_append_event[MAX_IN_FLIGHT] = {};

    // Pool bookkeeping, guarded by g_pool_mutex.
    bool in_use = false;
    uint64_t last_used_tick = 0;
    // Set when an exception mid-run leaves this instance's CUDA context/streams in an
    // unknown state; such instances must never be reused, only evicted and rebuilt.
    bool broken = false;

    bool matches(const GPUDecompressionConfig& c, int m, int dev) const {
        return gpu_slots == c.gpu_slots && reader_threads == c.reader_threads &&
               max_batch_size == c.max_batch_size &&
               max_residues == c.max_residues &&
               use_fused == c.use_fused_discretizer &&
               mode == m && device_id == dev;
    }
};

// The pool of resident pipeline instances (capped at the largest pool_size seen so
// far — a caller shrinking pool_size does not evict instances built under a larger
// one). Guarded by g_pool_mutex for structural changes (push/erase); individual
// instances' non-bookkeeping fields are only touched by whichever thread currently
// holds them checked out.
static std::vector<std::unique_ptr<PooledPipeline>> g_pool;

// Frees one instance's CUDA resources (streams/events/buffers), switching to its
// allocation device first if it differs from the caller's active device, then
// restores the caller's active device. Shared by eviction and full-pool release.
static void destroyPooledPipeline(PooledPipeline& p, int active_device) {
    const bool device_switch = p.device_id != active_device;
    if (device_switch) cudaSetDevice(p.device_id);
    CUDA_CHECK_FREE(cudaFreeHost(p.pdb_capture));
    for (cudaEvent_t& e : p.pdb_append_event) {
        if (e) { CUDA_CHECK(cudaEventDestroy(e)); e = nullptr; }
    }
    p.gpu_fc.clear();    // ~FoldcompGPUState destroys each slot's owned streams
    p.batch_pool.clear();
    if (device_switch) cudaSetDevice(active_device);
}

// Builds a brand-new instance in place (p is freshly default-constructed or just
// destroyPooledPipeline'd) for the given config/mode/device. Called with g_pool_mutex
// released (see buildPooledPipelineUnlocked) — the caller must already have set
// p.device_id and p.gpu_slots before releasing the lock, so this function leaves
// both untouched rather than re-writing them unsynchronized.
static void buildPooledPipeline(
    PooledPipeline& p, const GPUDecompressionConfig& config, int mode) {
    const int n_slots = config.gpu_slots;
    initGPUAminoAcidTable();
    const int preallocate_residues = config.max_batch_size * config.max_residues;
    const int preallocate_atoms = preallocate_residues * 10;
    p.gpu_fc.reserve(n_slots);
    for (int i = 0; i < n_slots; ++i) {
        auto fc = std::make_unique<Foldcomp>();
        fc->gpu.useFusedDisc = config.use_fused_discretizer;
        fc->gpu.initOwnedStreams();
        fc->gpu.gpu_disc_ctx.preallocate(config.max_batch_size, config.max_residues, /*max_sc_per_res=*/4,
            fc->gpu.streamA());
        fc->gpu.gpu_nerf_ctx.preallocate(config.max_batch_size, preallocate_residues, fc->gpu.streamA());
        // gpu_ac_buf (pinned AtomCoordinate host staging) is host-path only; the
        // pinned gpu_pdb_lines buffer grows lazily to the real batch size.
        const size_t ac_output_bytes =
            static_cast<size_t>(preallocate_atoms) * sizeof(AtomCoordinate);
        fc->gpu.gpu_sc_ctx.preallocate(config.max_batch_size, preallocate_residues, preallocate_atoms,
            ac_output_bytes, fc->gpu.streamA());
        if (mode == MODE_HOST) {
            fc->gpu.gpu_ac_buf.preallocate(
                static_cast<size_t>(preallocate_atoms) * sizeof(AtomCoordinate));
        }
        fc->gpu.gpu_bb_output_coords.reserve(static_cast<size_t>(preallocate_residues) * 9);
        p.gpu_fc.push_back(std::move(fc));
    }
    const int n_batches = n_slots + config.reader_threads;
    p.batch_pool.reserve(n_batches);
    for (int i = 0; i < n_batches; ++i) {
        auto batch = std::make_unique<Batch>();
        batch->data.max_batch_size = config.max_batch_size;
        batch->data.preallocate_pinned(
            static_cast<size_t>(config.max_batch_size) * config.max_residues,
            static_cast<size_t>(config.max_batch_size) * config.max_residues * 4);
        p.batch_pool.push_back(std::move(batch));
    }
    p.reader_threads = config.reader_threads;
    p.max_batch_size = config.max_batch_size;
    p.max_residues = config.max_residues;
    p.use_fused = config.use_fused_discretizer;
    p.mode = mode;
}

// Grows the shared CPU-sync event pool (gpu_stream.h) to cover the worst-case
// simultaneous in-flight count across every resident instance on active_device
// (sum of gpu_slots for that device only — the pool is keyed per-device).
// Must be called with g_pool_mutex held, after g_pool's membership changes.
static void growCpuSyncPoolForResidentInstances(int active_device) {
    int total_slots = 0;
    for (const std::unique_ptr<PooledPipeline>& p : g_pool)
        if (p->device_id == active_device) total_slots += p->gpu_slots;
    initGPUPipelineEvents(total_slots);
}

// Runs buildPooledPipeline on `p` with g_pool_mutex released, so the expensive part
// of a checkout (cudaMallocHost pinned preallocate, stream/event creation — all
// specific to active_device) never blocks a concurrent checkout/return on a
// different device. `p` must already be reserved (in_use = true, device_id and
// gpu_slots set to their final values) before the lock is released by the caller —
// those two fields are the only ones growCpuSyncPoolForResidentInstances reads for
// *every* pool entry regardless of in_use, so they must never be written while
// unlocked. On failure, the instance is marked broken (and freed for reuse via
// in_use = false) rather than left stuck in_use forever, then the exception is
// rethrown. Caller must hold the lock on entry and re-acquire it after return.
static void buildPooledPipelineUnlocked(
    std::unique_lock<std::mutex>& lock, PooledPipeline& p,
    const GPUDecompressionConfig& config, int mode) {
    lock.unlock();
    try {
        buildPooledPipeline(p, config, mode);
    } catch (...) {
        // Free whatever partial GPU/pinned-host resources buildPooledPipeline managed
        // to allocate before failing, right away rather than leaving them resident
        // until this broken instance is next evicted (which may be much later, or
        // never, if the pool never fills to capacity again).
        destroyPooledPipeline(p, p.device_id);
        lock.lock();
        p.broken = true;
        p.in_use = false;
        g_pool_cv.notify_all();
        throw;
    }
    lock.lock();
}

// Checks out a pipeline instance matching (config, mode, active_device), building or
// evicting-and-rebuilding one if needed, blocking if the pool is at capacity and
// every instance is busy. Returns with the instance marked in_use; caller must call
// returnPipeline() exactly once when done.
static PooledPipeline& checkoutPipeline(
    const GPUDecompressionConfig& config, int mode, int active_device) {
    std::unique_lock<std::mutex> lock(g_pool_mutex);
    while (true) {
        // 1. Idle instance with a matching fingerprint: reuse it.
        for (std::unique_ptr<PooledPipeline>& p : g_pool) {
            if (!p->in_use && !p->broken && p->matches(config, mode, active_device)) {
                p->in_use = true;
                return *p;
            }
        }
        // 2. Under capacity: build a new instance. Reserve it under the lock (so
        // concurrent callers see accurate pool size/idle state immediately), then
        // release the lock for the build itself — see buildPooledPipelineUnlocked.
        if (static_cast<int>(g_pool.size()) < config.pool_size) {
            auto p = std::make_unique<PooledPipeline>();
            p->in_use = true;
            p->device_id = active_device;
            p->gpu_slots = config.gpu_slots;
            PooledPipeline& ref = *p;
            g_pool.push_back(std::move(p));
            buildPooledPipelineUnlocked(lock, ref, config, mode);
            // Must run after the build above: it sums gpu_slots across g_pool,
            // which needs to already include this new instance.
            growCpuSyncPoolForResidentInstances(active_device);
            return ref;
        }
        // 3. At capacity: evict the least-recently-used idle instance, if any.
        PooledPipeline* victim = nullptr;
        for (std::unique_ptr<PooledPipeline>& p : g_pool) {
            if (p->in_use) continue;
            if (p->broken) { victim = p.get(); break; } // broken instances are unusable; evict on sight
            if (!victim || p->last_used_tick < victim->last_used_tick) victim = p.get();
        }
        if (victim) {
            destroyPooledPipeline(*victim, active_device);
            *victim = PooledPipeline(); // reset to a fresh default state before rebuild
            victim->in_use = true;
            victim->device_id = active_device;
            victim->gpu_slots = config.gpu_slots;
            buildPooledPipelineUnlocked(lock, *victim, config, mode);
            growCpuSyncPoolForResidentInstances(active_device);
            return *victim;
        }
        // 4. Every instance busy: wait for one to be returned, then retry. Bounded so a
        // stuck/hung in-flight run (e.g. a wedged CUDA stream) fails callers fast instead
        // of blocking them forever.
        constexpr std::chrono::seconds kCheckoutTimeout{5};
        if (g_pool_cv.wait_for(lock, kCheckoutTimeout) == std::cv_status::timeout) {
            throw std::runtime_error(
                "GPU decompression pipeline pool exhausted: timed out after 5s waiting "
                "for a free slot (all pool_size instances busy or stuck)");
        }
    }
}

static void returnPipeline(PooledPipeline& p) {
    std::lock_guard<std::mutex> lock(g_pool_mutex);
    p.in_use = false;
    p.last_used_tick = ++g_pool_tick;
    g_pool_cv.notify_all();
}

GPUPipelineLease::~GPUPipelineLease() {
    if (pipe_) returnPipeline(*static_cast<PooledPipeline*>(pipe_));
}

GPUPipelineLease& GPUPipelineLease::operator=(GPUPipelineLease&& other) noexcept {
    if (this != &other) {
        if (pipe_) returnPipeline(*static_cast<PooledPipeline*>(pipe_));
        pipe_ = other.pipe_;
        other.pipe_ = nullptr;
    }
    return *this;
}

const char* GPUPipelineLease::pdbCaptureBase() const {
    return pipe_ ? static_cast<PooledPipeline*>(pipe_)->pdb_capture : nullptr;
}

static void recordAppendEvent(cudaEvent_t& slot_event, cudaStream_t stream) {
    if (!slot_event)
        CUDA_CHECK(cudaEventCreateWithFlags(&slot_event, cudaEventDisableTiming));
    CUDA_CHECK(cudaEventRecord(slot_event, stream));
}

// Wait for every recorded per-slot append event instead of the whole device.
static void waitAppendEvents(cudaEvent_t (&events)[MAX_IN_FLIGHT]) {
    for (cudaEvent_t e : events)
        if (e) CUDA_CHECK(cudaEventSynchronize(e));
}

// Ensure pipe's PDB capture buffer holds at least `need_atoms` lines, preserving
// already appended data on growth. Growth is rare (buffer stabilizes at the
// high-water mark across calls that reuse this instance). On growth, waits only on
// the per-slot append events (see pdb_append_event above) instead of the whole
// device, so a resize doesn't stall slots that don't write into this buffer.
static char* ensurePdbCapture(PooledPipeline& pipe, size_t need_atoms) {
    const size_t need_bytes = need_atoms * PDB_ATOM_LINE_LEN;
    if (need_bytes <= pipe.pdb_capacity) return pipe.pdb_capture;
    const size_t new_cap = std::max(need_bytes, pipe.pdb_capacity * 2);
    char* np = nullptr;
    if (cudaMallocHost(&np, new_cap) != cudaSuccess)
        throw std::runtime_error("pdb line capture cudaMallocHost failed");
    if (pipe.pdb_capture && pipe.pdb_used > 0) {
        waitAppendEvents(pipe.pdb_append_event); // let pending D→H land in the old buffer first
        std::memcpy(np, pipe.pdb_capture, pipe.pdb_used * PDB_ATOM_LINE_LEN);
    }
    CUDA_CHECK_FREE(cudaFreeHost(pipe.pdb_capture));
    pipe.pdb_capture = np;
    pipe.pdb_capacity = new_cap;
    return np;
}

void releaseGPUPipelineCache() {
    std::unique_lock<std::mutex> lock(g_pool_mutex);
    // Bounded like checkoutPipeline's wait above: a stuck/hung in-flight run must not
    // hang this call (and thus the calling Python thread) forever.
    constexpr std::chrono::seconds kReleaseTimeout{5};
    auto all_idle = [] {
        for (const std::unique_ptr<PooledPipeline>& p : g_pool)
            if (p->in_use) return false;
        return true;
    };
    if (!g_pool_cv.wait_for(lock, kReleaseTimeout, all_idle)) {
        throw std::runtime_error(
            "GPU decompression pipeline release timed out after 5s: a pooled instance "
            "is still marked in-use (busy or stuck)");
    }
    int active_device = 0;
    cudaGetDevice(&active_device);
    for (std::unique_ptr<PooledPipeline>& p : g_pool) destroyPooledPipeline(*p, active_device);
    g_pool.clear();
    g_pool_tick = 0;
    destroyGPUPipelineEvents();
}

bool GPUDecompressionConfig::validate(std::string& error) const {
    if (gpu_slots < 1) {
        error = "GPU pipeline gpu_slots must be positive (got " +
            std::to_string(gpu_slots) + ")";
        return false;
    }
    if (gpu_slots > MAX_IN_FLIGHT) {
        error = "GPU pipeline gpu_slots exceeds MAX_IN_FLIGHT (" +
            std::to_string(MAX_IN_FLIGHT) + ")";
        return false;
    }
    if (reader_threads < 1) {
        error = "GPU pipeline reader_threads must be positive (got " +
            std::to_string(reader_threads) + ")";
        return false;
    }
    if (max_batch_size < 1) {
        error = "GPU pipeline max_batch_size must be positive (got " +
            std::to_string(max_batch_size) + ")";
        return false;
    }
    if (max_residues < 1) {
        error = "GPU pipeline max_residues must be positive (got " +
            std::to_string(max_residues) + ")";
        return false;
    }
    if (pool_size < 1) {
        error = "GPU pipeline pool_size must be positive (got " +
            std::to_string(pool_size) + ")";
        return false;
    }
    if (max_batch_size > FoldcompGPUState::MAX_STRUCTS_PER_BATCH) {
        error = "GPU pipeline max_batch_size exceeds MAX_STRUCTS_PER_BATCH (" +
            std::to_string(FoldcompGPUState::MAX_STRUCTS_PER_BATCH) + ")";
        return false;
    }
    return true;
}

static bool runGPUDecompressionPipelineImpl(
    Processor& processor,
    const GPUDecompressionConfig& config,
    const GPUDecompressionOutputFn& output,
    std::string& error,
    GPUPipelineLease* out_lease) {
    const bool pdb_mode = config.emit_device_pdb;
    const int mode = pdb_mode ? MODE_PDB : MODE_HOST;

    error.clear();
    if (!config.validate(error)) return false;
    if (!cudaRuntimeAvailable(&error)) {
        error = "CUDA runtime is unavailable: " + error;
        return false;
    }

    // Check out a pipeline instance for this config/mode/device; only the per-call
    // queues, reader thread, and parsers below are transient. The instance is
    // returned to the pool when `guard` is destroyed, unless ownership is handed to
    // `out_lease` at the end of this function (see below).
    int active_device = 0;
    cudaGetDevice(&active_device);
    PooledPipeline& pipe = checkoutPipeline(config, mode, active_device);
    struct ReturnGuard {
        PooledPipeline* p;
        bool leased = false;
        ~ReturnGuard() { if (p && !leased) returnPipeline(*p); }
    } guard{&pipe};
    const int n_slots = pipe.gpu_slots;
    std::vector<std::unique_ptr<Foldcomp>>& gpu_fc = pipe.gpu_fc;
    std::vector<std::unique_ptr<Batch>>& batch_pool = pipe.batch_pool;

    std::queue<Batch*> free_q;
    std::queue<Batch*> ready_q;
    std::queue<std::shared_ptr<LogicalAssembly>> completed_q;
    std::mutex free_mtx, ready_mtx, error_mtx;
    std::condition_variable cv_free, cv_ready;
    // Reset the persisted batches to a clean state and stage them as free.
    for (std::unique_ptr<Batch>& b : batch_pool) {
        b->data.clear();
        b->targets.clear();
        free_q.push(b.get());
    }
    if (pdb_mode) pipe.pdb_used = 0; // start a fresh line capture

    std::vector<Foldcomp> parsers(config.reader_threads);
    std::vector<Batch*> thread_batches(config.reader_threads, nullptr);
    std::atomic<bool> reader_done{false};
    std::atomic<bool> read_failed{false};
    std::atomic<bool> cancelled{false};

    auto set_error = [&](const std::string& message) {
        read_failed.store(true, std::memory_order_release);
        std::lock_guard<std::mutex> lock(error_mtx);
        if (error.empty()) error = message;
    };

    auto acquire_batch = [&]() -> Batch* {
        std::unique_lock<std::mutex> lock(free_mtx);
        cv_free.wait(lock, [&] {
            return !free_q.empty() || cancelled.load(std::memory_order_acquire);
        });
        if (cancelled.load(std::memory_order_acquire)) return nullptr;
        Batch* batch = free_q.front();
        free_q.pop();
        return batch;
    };

    auto publish_batch = [&](Batch*& batch) {
        if (!batch) return;
        if (batch->data.current_batch_count == 0) {
            std::lock_guard<std::mutex> lock(free_mtx);
            free_q.push(batch);
            cv_free.notify_one();
            batch = nullptr;
            return;
        }
        batch->data.finalize();
        {
            std::lock_guard<std::mutex> lock(ready_mtx);
            ready_q.push(batch);
        }
        batch = nullptr;
        cv_ready.notify_one();
    };

    process_entry_func reader = [&](const char* name, const char* data, size_t size) -> bool {
        if (!name) return true;
        if (cancelled.load(std::memory_order_acquire)) return false;
        try {
            int tid = omp_get_thread_num();
            if (tid < 0 || tid >= config.reader_threads) tid = 0;
            ParsedInput parsed;
            std::string parse_error;
            if (!parseInput(parsers[tid], name, data, size, config, parsed, parse_error)) {
                set_error(std::string(name) + ": " + parse_error);
                return false;
            }

            if (parsed.gpu_fragments.empty()) {
                std::lock_guard<std::mutex> lock(ready_mtx);
                completed_q.push(std::move(parsed.assembly));
                cv_ready.notify_one();
                return true;
            }

            Batch*& batch = thread_batches[tid];
            for (ParsedGPUFragment& fragment : parsed.gpu_fragments) {
                if (!batch) batch = acquire_batch();
                if (!batch) return false;
                batch->data.append_blob(std::move(fragment.blob));
                batch->targets.push_back(std::move(fragment.target));
                if (batch->data.is_full()) publish_batch(batch);
            }
            return true;
        } catch (const std::exception& exception) {
            set_error(std::string(name) + ": " + exception.what());
            return false;
        } catch (...) {
            set_error(std::string(name) + ": unknown reader failure");
            return false;
        }
    };

    std::thread reader_thread([&] {
        try {
            processor.run(reader, config.reader_threads);
            if (!cancelled.load(std::memory_order_acquire)) {
                for (Batch*& batch : thread_batches) publish_batch(batch);
            }
        } catch (const std::exception& exception) {
            set_error(exception.what());
        } catch (...) {
            set_error("unknown input processor failure");
        }
        reader_done.store(true, std::memory_order_release);
        cv_ready.notify_all();
    });

    struct Pending {
        Batch* batch = nullptr;
        int fc_index = 0;
        int count = 0;
        std::vector<int> atom_counts;
        AtomCoordinate* atoms = nullptr;
        cudaEvent_t event = nullptr;
        size_t pdb_base = 0; // atom offset of this sub-batch in pipe.pdb_capture
    };
    std::deque<Pending> pending;
    int batch_index = 0;
    bool output_ok = true;

    auto emit_host = [&](const std::shared_ptr<LogicalAssembly>& assembly) {
        if (!output_ok || assembly->emitted) return;
        assembly->emitted = true;
        GPUDecompressionResult result;
        result.name = assembly->name;
        result.title = assembly->title;
        result.segments = std::move(assembly->slots);
        if (!output(std::move(result))) {
            output_ok = false;
            set_error(assembly->name + ": output callback failed");
        }
    };

    auto ensure_slot_pdb = [](LogicalAssembly& a) {
        if (a.slot_pdb_off.empty()) {
            a.slot_pdb_off.assign(a.slots.size(), SIZE_MAX);
            a.slot_pdb_meta.assign(a.slots.size(), PdbLineSegment{0, ChainId(), 0});
        }
    };

    // emit_device_pdb: concatenate an assembly's per-slot GPU-formatted ATOM lines
    // (formatting raw/CPU-fallback fragments on the host) plus per-segment metadata,
    // delivered through the host output callback.
    auto emit_pdb_lines = [&](const std::shared_ptr<LogicalAssembly>& assembly) {
        if (!output_ok || assembly->emitted) return;
        assembly->emitted = true;
        ensure_slot_pdb(*assembly);
        GPUDecompressionResult result;
        result.name = assembly->name;
        result.title = assembly->title;
        result.line_segments.reserve(assembly->slots.size());
        result.pdb_frags.reserve(assembly->slots.size());
        for (size_t i = 0; i < assembly->slots.size(); ++i) {
            if (assembly->slot_pdb_off[i] != SIZE_MAX) {
                result.line_segments.push_back(assembly->slot_pdb_meta[i]);
                result.pdb_frags.push_back(
                    {assembly->slot_pdb_off[i], assembly->slot_pdb_meta[i].n_atoms});
            } else {
                // Raw/CPU-fallback fragment: format its atoms into the shared line
                // capture buffer (host, rare) and record its offset.
                const std::vector<AtomCoordinate>& atoms = assembly->slots[i].atoms;
                const int n = static_cast<int>(atoms.size());
                ChainId chain = (n == 0 || atoms[0].chain.empty()) ? ChainId(" ") : atoms[0].chain;
                const size_t off = pipe.pdb_used;
                char* cap = ensurePdbCapture(pipe, off + static_cast<size_t>(n));
                for (int k = 0; k < n; ++k) {
                    formatAtomLineHost(atoms[k],
                        cap + (off + static_cast<size_t>(k)) * PDB_ATOM_LINE_LEN);
                }
                pipe.pdb_used += static_cast<size_t>(n);
                result.line_segments.push_back(
                    PdbLineSegment{assembly->slots[i].model, chain, n});
                result.pdb_frags.push_back({off, n});
            }
        }
        if (!output(std::move(result))) {
            output_ok = false;
            set_error(assembly->name + ": output callback failed");
        }
    };

    auto emit = [&](const std::shared_ptr<LogicalAssembly>& assembly) {
        if (pdb_mode) emit_pdb_lines(assembly);
        else emit_host(assembly);
    };

    auto drain = [&] {
        Pending current = std::move(pending.front());
        pending.pop_front();
        CUDA_CHECK(cudaEventSynchronize(current.event));
        releaseCPUSyncEvent(current.event);
        if (pdb_mode) {
            Foldcomp* fc = gpu_fc[current.fc_index].get();
            // Lines were D→H'd directly into this instance's capture buffer at pdb_base.
            char* sub_base = pipe.pdb_capture +
                current.pdb_base * PDB_ATOM_LINE_LEN;
            fc->finish_pending_oxt_lines(sub_base); // patch OXT slots in place
            int atom_offset = 0;
            for (int i = 0; i < current.count; ++i) {
                const int n_atoms = current.atom_counts[i];
                FragmentTarget& target = current.batch->targets[i];
                LogicalAssembly& a = *target.assembly;
                ensure_slot_pdb(a);
                if (output_ok && n_atoms > 0) {
                    const size_t frag_off = current.pdb_base + atom_offset;
                    char* seg = pipe.pdb_capture + frag_off * PDB_ATOM_LINE_LEN;
                    ChainId chain;
                    if (target.override_chain && !target.chain.empty()) {
                        // Full-width chain id (may be >1 char, e.g. mmCIF auth_asym_id)
                        // for the boundary check below; the PDB format itself only has
                        // one column for chain, so that column still gets just chain[0].
                        chain = target.chain;
                        char col21 = target.chain[0];
                        for (int k = 0; k < n_atoms; ++k) {
                            seg[static_cast<size_t>(k) * PDB_ATOM_LINE_LEN + 21] = col21;
                        }
                    } else {
                        // Non-container FCZ fragments carry a single-char chain at the
                        // source (foldcomp.h: char chain), so the formatted line's
                        // column-21 byte already has the full chain identity.
                        chain = ChainId(std::string(1, seg[21]));
                    }
                    // Record the fragment's offset into the shared buffer (no copy).
                    a.slot_pdb_off[target.slot] = frag_off;
                    a.slot_pdb_meta[target.slot] =
                        PdbLineSegment{target.model, chain, n_atoms};
                }
                atom_offset += n_atoms;
                a.pending_gpu--;
                if (a.pending_gpu == 0) emit(target.assembly);
            }
        } else {
            gpu_fc[current.fc_index]->finish_pending_oxt(current.atoms);
            int atom_offset = 0;
            for (int i = 0; i < current.count; ++i) {
                const int n_atoms = current.atom_counts[i];
                FragmentTarget& target = current.batch->targets[i];
                std::vector<AtomCoordinate> atoms(
                    current.atoms + atom_offset,
                    current.atoms + atom_offset + n_atoms);
                atom_offset += n_atoms;
                for (AtomCoordinate& atom : atoms) {
                    atom.model = target.model;
                    if (target.override_chain) atom.chain = target.chain;
                }
                target.assembly->slots[target.slot] = {target.model, std::move(atoms)};
                target.assembly->pending_gpu--;
                if (target.assembly->pending_gpu == 0) emit(target.assembly);
            }
        }
        current.batch->data.clear();
        current.batch->targets.clear();
        {
            std::lock_guard<std::mutex> lock(free_mtx);
            free_q.push(current.batch);
        }
        cv_free.notify_one();
    };

    try {
        while (true) {
            Batch* batch = nullptr;
            std::shared_ptr<LogicalAssembly> completed;
            {
                std::unique_lock<std::mutex> lock(ready_mtx);
                cv_ready.wait(lock, [&] {
                    return !ready_q.empty() || !completed_q.empty() ||
                           reader_done.load(std::memory_order_acquire);
                });
                if (!completed_q.empty()) {
                    completed = std::move(completed_q.front());
                    completed_q.pop();
                } else if (!ready_q.empty()) {
                    batch = ready_q.front();
                    ready_q.pop();
                } else {
                    lock.unlock();
                    while (!pending.empty()) drain();
                    break;
                }
            }
            if (completed) {
                emit(completed);
                continue;
            }

            const int fc_index = batch_index++ % n_slots;
            Pending item;
            item.batch = batch;
            item.fc_index = fc_index;
            item.count = batch->data.current_batch_count;
            // Placeholder until decompress_batch_async below fills in the real
            // per-structure counts; only used as-is if current_batch_count == 0
            // (empty batch, decompress_batch_async returns before touching it).
            item.atom_counts = batch->data.nAtom_batch;
            const int total_atoms = batch->data.cached_total_atoms;
            item.atoms = (mode == MODE_HOST)
                ? gpu_fc[fc_index]->gpu.gpu_ac_buf.ensure_n<AtomCoordinate>(total_atoms)
                : nullptr;
            // PDB mode: reserve this sub-batch's slice of the shared line buffer and
            // have the async D→H land the formatted lines directly there.
            char* pdb_out = nullptr;
            if (pdb_mode) {
                item.pdb_base = pipe.pdb_used;
                char* cap = ensurePdbCapture(pipe, pipe.pdb_used + total_atoms);
                pdb_out = cap + pipe.pdb_used * PDB_ATOM_LINE_LEN;
                pipe.pdb_used += total_atoms;
            }
            item.event = acquireCPUSyncEvent();
            AtomCoordinate* const item_atoms = item.atoms;
            cudaEvent_t const item_event = item.event;
            std::swap(gpu_fc[fc_index]->gpu.batchData, batch->data);
            // Push into `pending` before the async dispatch below: decompress_batch_async
            // and recordAppendEvent can both throw via CUDA_CHECK, and the catch blocks'
            // drain loop only reclaims events from `pending` — pushing first ensures
            // item.event (already acquired above) is always reachable for
            // releaseCPUSyncEvent instead of leaking one pool event per such failure.
            pending.push_back(std::move(item));
            // Fill pending.back().atom_counts in place with the real per-structure
            // counts decompress_batch_async computes internally (may differ from the
            // nAtom_batch placeholder above — see the out_atom_counts doc comment on
            // decompress_batch_async). pending.back() stays valid: nothing else grows
            // `pending` between this push_back and the call below.
            gpu_fc[fc_index]->decompress_batch_async(
                item_atoms, total_atoms, item_event, /*device_emit=*/false, pdb_mode, pdb_out,
                &pending.back().atom_counts);
            if (pdb_mode) {
                // decompress_batch_async's internal D→H write into pdb_out lands on
                // this same stream; record after it's enqueued so a growth triggered
                // by a later dispatch on another slot waits for it precisely.
                recordAppendEvent(pipe.pdb_append_event[fc_index], gpu_fc[fc_index]->gpu.streamA());
            }
            if (static_cast<int>(pending.size()) >= n_slots) drain();
        }
    } catch (const std::exception& exception) {
        set_error(exception.what());
        cancelled.store(true, std::memory_order_release);
        cv_free.notify_all();
        cv_ready.notify_all();
        if (cudaDeviceSynchronize() != cudaSuccess) pipe.broken = true;
        while (!pending.empty()) {
            cudaEventSynchronize(pending.front().event);
            releaseCPUSyncEvent(pending.front().event);
            pending.pop_front();
        }
    } catch (...) {
        set_error("unknown GPU pipeline failure");
        cancelled.store(true, std::memory_order_release);
        cv_free.notify_all();
        cv_ready.notify_all();
        if (cudaDeviceSynchronize() != cudaSuccess) pipe.broken = true;
        while (!pending.empty()) {
            cudaEventSynchronize(pending.front().event);
            releaseCPUSyncEvent(pending.front().event);
            pending.pop_front();
        }
    }

    reader_thread.join();
    // Ensure all async PDB line D→H writes have landed before the caller reads this
    // instance's capture buffer. This is also the only place a deferred error from a
    // kernel launched during the pipeline (sidechain/backbone fault, out-of-bounds
    // device write, etc.) surfaces on the otherwise-successful path, so it needs the
    // same broken/output_ok handling as the catch blocks above — otherwise the caller
    // could read a silently corrupt capture buffer.
    if (pdb_mode) {
        if (cudaDeviceSynchronize() != cudaSuccess) {
            pipe.broken = true;
            output_ok = false;
            set_error("cudaDeviceSynchronize failed after GPU pipeline run");
        }
    }

    // Hand ownership of the instance to the caller's lease instead of returning it
    // to the pool, so the caller can safely read pipe.pdb_capture (via
    // GPUPipelineLease) after this function returns.
    if (out_lease) {
        out_lease->pipe_ = &pipe;
        guard.leased = true;
    }
    return !read_failed.load(std::memory_order_acquire) && output_ok;
}

bool runGPUDecompressionPipeline(
    Processor& processor,
    const GPUDecompressionConfig& config,
    const GPUDecompressionOutputFn& output,
    std::string& error,
    GPUPipelineLease* out_lease) {
    return runGPUDecompressionPipelineImpl(processor, config, output, error, out_lease);
}

#endif
