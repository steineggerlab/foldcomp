// SPDX-License-Identifier: MIT
//
// C++ counterpart to test/benchmark_cpu_gpu.py. Measures, for batches of
// 8/16/32/64/128/256 identical structures, the three core Foldcomp decompression
// paths — calling the SAME internal functions the Python binding uses, but
// without the Python/PyUnicode/GIL layer:
//
//   * CPU      - decodeStructureToPDB per structure (CPU decompress -> PDB text)
//   * GPUhost  - runGPUDecompressionPipeline (host format) + writeSegmentsToPDB
//   * GPUfmt   - pipeline with device PDB formatting + parallel line stitch
//
// Nothing is written to disk. Methodology mirrors the Python harness: median
// wall time, one warmup, warmed GPU pipeline cache (first-call setup amortized),
// and per-iteration NVTX ranges with the same labels ("GPUfmt N=32/iter3", ...)
// so a C++ trace lines up directly against a Python trace.
//
// Build:  configure with -DBUILD_CUDA=ON (executable build) -> `foldcomp_bench`.
// Run:    ./foldcomp_bench [--iters N] [--sizes 8,32,...] [--structure test/test.pdb]

#include <cuda_profiler_api.h>
#include <cuda_runtime.h>
#include <nvtx3/nvToolsExt.h>
#include <nvtx3/nvtx3.hpp>

#include "atom_coordinate.h"          // PDB_ATOM_LINE_LEN
#include "gpu_decompression_pipeline.h"
#include "input_processor.h"
#include "structure_codec.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

namespace {

using Clock = std::chrono::steady_clock;

// Foldcomp GPU-path defaults (match the Python binding's decompress_batch* args).
int g_max_batch_size = 100;
int g_max_residues = 3000;

std::string readFile(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) {
        fprintf(stderr, "[error] cannot open %s\n", path.c_str());
        std::exit(1);
    }
    return std::string((std::istreambuf_iterator<char>(f)),
                       std::istreambuf_iterator<char>());
}

// Compress a PDB/CIF file to an FCZ container, or pass an FCZ file through
// unchanged — the analog of Python's foldcomp.compress(test.pdb).
std::string toFcz(const std::string& path) {
    std::string raw = readFile(path);
    if (hasContainerMagic(raw.data(), raw.size())) {
        return raw;  // already an FCZ container
    }
    const char* fmt = (path.size() >= 4 &&
                       (path.compare(path.size() - 4, 4, ".cif") == 0 ||
                        path.compare(path.size() - 4, 4, ".CIF") == 0))
                          ? "cif" : "pdb";
    std::string title = path;
    const size_t slash = title.find_last_of('/');
    if (slash != std::string::npos) title = title.substr(slash + 1);
    std::string out;
    int rc = encodeStructureToFoldcomp(
        title, raw.data(), raw.size(), fmt,
        /*anchorResidueThreshold=*/25,
        /*maxBackboneRmsd=*/std::numeric_limits<float>::infinity(), out);
    if (rc != 0) {
        fprintf(stderr, "[error] failed to compress %s (rc=%d)\n", path.c_str(), rc);
        std::exit(1);
    }
    return out;
}

double medianMs(std::vector<double>& v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

// Median wall time (ms) of `fn` over `iters` runs after `warmup` unmeasured
// runs. Each run is wrapped in an NVTX range labelled "<label>/iterK" (or
// "/warmupW") to match the Python harness.
template <class F>
double timeit(F&& fn, int iters, const std::string& label, int warmup = 1) {
    for (int w = 0; w < warmup; ++w) {
        nvtx3::scoped_range r{label + "/warmup" + std::to_string(w)};
        fn();
    }
    std::vector<double> samples;
    samples.reserve(static_cast<size_t>(iters));
    for (int k = 0; k < iters; ++k) {
        nvtx3::scoped_range r{label + "/iter" + std::to_string(k)};
        const auto t0 = Clock::now();
        fn();
        const auto t1 = Clock::now();
        samples.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
    }
    return medianMs(samples);
}

std::vector<MemoryProcessor::Entry> makeEntries(const std::vector<std::string>& batch) {
    std::vector<MemoryProcessor::Entry> entries;
    entries.reserve(batch.size());
    for (size_t i = 0; i < batch.size(); ++i) {
        entries.emplace_back(std::to_string(i),
                             std::vector<char>(batch[i].begin(), batch[i].end()));
    }
    return entries;
}

GPUDecompressionConfig baseConfig() {
    GPUDecompressionConfig cfg;
    cfg.reader_threads = 1;  // matches the Python binding
    cfg.max_batch_size = g_max_batch_size;
    cfg.max_residues = g_max_residues;
    return cfg;
}

// --- CPU: decodeStructureToPDB per structure (CPU decompress -> PDB text) -----
void cpuDecompress(const std::vector<std::string>& batch,
                   std::vector<std::string>& out) {
    out.assign(batch.size(), std::string());
    for (size_t i = 0; i < batch.size(); ++i) {
        std::string title;
        decodeStructureToPDB(batch[i].data(), batch[i].size(),
                             /*useAltOrder=*/false, title, out[i]);
    }
}

// --- GPUhost: pipeline (host format) + writeSegmentsToPDB per structure -------
bool gpuHostDecompress(const std::vector<std::string>& batch,
                       std::vector<std::string>& out) {
    auto entries = makeEntries(batch);
    MemoryProcessor processor(entries);
    GPUDecompressionConfig cfg = baseConfig();
    out.assign(batch.size(), std::string());
    auto collect = [&](GPUDecompressionResult&& r) -> bool {
        const size_t idx = static_cast<size_t>(std::stoul(r.name));
        if (idx >= out.size()) return false;
        writeSegmentsToPDB(r.segments, r.title, out[idx]);
        return true;
    };
    std::string error;
    if (!runGPUDecompressionPipeline(processor, cfg, collect, error)) {
        fprintf(stderr, "[error] GPUhost: %s\n", error.c_str());
        return false;
    }
    return true;
}

// --- GPUfmt: pipeline with device PDB formatting + parallel line stitch -------
//     (replicates foldcomp.cxx's deferred parallel stitch)
bool gpuFmtDecompress(const std::vector<std::string>& batch,
                      std::vector<std::string>& out) {
    auto entries = makeEntries(batch);
    MemoryProcessor processor(entries);
    GPUDecompressionConfig cfg = baseConfig();
    cfg.emit_device_pdb = true;
    const size_t n = batch.size();
    out.assign(n, std::string());
    std::vector<GPUDecompressionResult> staging(n);
    auto collect = [&](GPUDecompressionResult&& r) -> bool {
        const size_t idx = static_cast<size_t>(std::stoul(r.name));
        if (idx >= staging.size()) return false;
        staging[idx] = std::move(r);
        return true;
    };
    std::string error;
    GPUPipelineLease lease;
    {
        nvtx3::scoped_range rr{"GPU pipeline (run)"};
        if (!runGPUDecompressionPipeline(processor, cfg, collect, error, &lease)) {
            fprintf(stderr, "[error] GPUfmt: %s\n", error.c_str());
            return false;
        }
    }
    {
        nvtx3::scoped_range rr{"GPU PDB stitch (host)"};
        const char* base = lease.pdbCaptureBase();
        const long count = static_cast<long>(n);
#pragma omp parallel for schedule(dynamic)
        for (long i = 0; i < count; ++i) {
            const GPUDecompressionResult& r = staging[static_cast<size_t>(i)];
            std::string& o = out[static_cast<size_t>(i)];
            if (r.pdb_frags.size() == 1) {
                writeLineSegmentsToPDB(
                    r.line_segments, r.title,
                    base + r.pdb_frags[0].offset * PDB_ATOM_LINE_LEN, o);
            } else {
                size_t tot = 0;
                for (const GPUCoordFragment& f : r.pdb_frags)
                    tot += static_cast<size_t>(f.n_atoms);
                std::string tmp(tot * PDB_ATOM_LINE_LEN, '\0');
                size_t off = 0;
                for (const GPUCoordFragment& f : r.pdb_frags) {
                    std::memcpy(&tmp[off * PDB_ATOM_LINE_LEN],
                                base + f.offset * PDB_ATOM_LINE_LEN,
                                static_cast<size_t>(f.n_atoms) * PDB_ATOM_LINE_LEN);
                    off += static_cast<size_t>(f.n_atoms);
                }
                writeLineSegmentsToPDB(r.line_segments, r.title, tmp.data(), o);
            }
        }
    }
    return true;
}

// Lightweight correctness check (analog of the Python verify_single): the GPU
// device-format text must be byte-identical to the GPU host-format text.
void verifySingle(const std::string& base_fcz) {
    std::vector<std::string> host, fmt;
    if (!gpuHostDecompress({base_fcz}, host)) return;
    if (!gpuFmtDecompress({base_fcz}, fmt)) return;
    std::vector<std::string> cpu;
    cpuDecompress({base_fcz}, cpu);
    const bool identical = (host.size() == 1 && fmt.size() == 1 && host[0] == fmt[0]);
    const size_t atoms = cpu.empty() ? 0 : (cpu[0].size() / (PDB_ATOM_LINE_LEN));
    printf("Single-file check: GPUfmt == GPUhost text: %s  (CPU text bytes=%zu)\n",
           identical ? "OK" : "FAIL", cpu.empty() ? size_t(0) : cpu[0].size());
    (void)atoms;
}

struct Row {
    int n_struct;
    size_t n_atoms;
    double cpu, gpu_host, gpu_fmt;  // ms
};

Row benchBatch(const std::vector<std::string>& batch, int iters) {
    const int n = static_cast<int>(batch.size());
    std::vector<std::string> textOut;

    const double t_cpu = timeit(
        [&] { cpuDecompress(batch, textOut); }, iters, "CPU N=" + std::to_string(n));
    size_t n_atoms = 0;
    for (const std::string& text : textOut) n_atoms += text.size() / PDB_ATOM_LINE_LEN;

    const double t_gpu_host = timeit(
        [&] { gpuHostDecompress(batch, textOut); }, iters, "GPUhost N=" + std::to_string(n));
    const double t_gpu_fmt = timeit(
        [&] { gpuFmtDecompress(batch, textOut); }, iters, "GPUfmt N=" + std::to_string(n));

    return Row{n, n_atoms, t_cpu, t_gpu_host, t_gpu_fmt};
}

void printRow(const Row& r) {
    const double speedup = (r.gpu_fmt > 0) ? r.cpu / r.gpu_fmt : 0.0;
    const double rate = (r.gpu_fmt > 0) ? r.n_struct / (r.gpu_fmt / 1e3) : 0.0;
    printf("%5d | %8zu | %9.2f | %10.2f | %9.2f | %6.2fx | %9.0f\n",
           r.n_struct, r.n_atoms, r.cpu, r.gpu_host, r.gpu_fmt, speedup, rate);
}

std::vector<int> parseSizes(const std::string& s) {
    std::vector<int> out;
    size_t i = 0;
    while (i < s.size()) {
        size_t j = s.find(',', i);
        if (j == std::string::npos) j = s.size();
        const std::string tok = s.substr(i, j - i);
        if (!tok.empty()) {
            int v = std::stoi(tok);
            if (v <= 0) {
                fprintf(stderr, "[error] --sizes entries must be positive, got %d\n", v);
                std::exit(2);
            }
            out.push_back(v);
        }
        i = j + 1;
    }
    return out;
}

int parsePositiveInt(const char* name, const std::string& value) {
    int v = std::stoi(value);
    if (v <= 0) {
        fprintf(stderr, "[error] %s must be positive, got %d\n", name, v);
        std::exit(2);
    }
    return v;
}

int parseNonNegativeInt(const char* name, const std::string& value) {
    int v = std::stoi(value);
    if (v < 0) {
        fprintf(stderr, "[error] %s must be non-negative, got %d\n", name, v);
        std::exit(2);
    }
    return v;
}

// Real-data bench mode: loads every file under `inputDir` into RAM once, then
// runs `warmup + samples` full-pipeline passes (CPU decodeStructureToPDB or
// GPU host-format pipeline) over the whole preloaded set, discarding output.
// Each measured pass is bracketed by cudaProfilerStart/Stop and an NVTX range
// so it can be captured with `nsys profile --capture-range=cudaProfilerApi`.
// Emits one JSON line with per-sample wall times, mirroring the retired
// `foldcomp decompress --bench-warmup/--bench-samples` CLI mode.
void runInputDirBench(const std::string& inputDir, bool recursive, bool cpuMode,
                      int warmup, int samples) {
    DirectoryProcessor proc(inputDir, recursive, /*use_mmap=*/true, /*sort_by_length=*/false);
    std::vector<MemoryProcessor::Entry> preloaded;
    proc.run([&](const char* name, const char* buf, size_t sz) -> bool {
        if (!name) return true;
        preloaded.push_back({std::string(name), std::vector<char>(buf, buf + sz)});
        return true;
    }, 1);
    if (preloaded.empty()) {
        fprintf(stderr, "[error] no files found under %s\n", inputDir.c_str());
        std::exit(1);
    }

    size_t batch_bytes = 0;
    for (const auto& e : preloaded) batch_bytes += e.second.size();
    fprintf(stderr, "[bench] pre-loaded %zu files (%.2f MB)\n",
            preloaded.size(), batch_bytes / 1e6);

    std::vector<std::string> batch;
    batch.reserve(preloaded.size());
    for (const auto& e : preloaded) batch.emplace_back(e.second.begin(), e.second.end());

    const int total_iters = warmup + samples;
    std::vector<double> times;
    times.reserve(static_cast<size_t>(samples));
    std::vector<std::string> textOut;

    for (int iter = 0; iter < total_iters; ++iter) {
        const bool measured = iter >= warmup;
        nvtxRangeId_t rng = 0;
        if (measured) {
            cudaProfilerStart();
            rng = nvtxRangeStartA(("sample_" + std::to_string(iter - warmup) +
                                   (cpuMode ? "_cpu" : "_gpu")).c_str());
        }
        const auto t0 = Clock::now();
        if (cpuMode) {
            cpuDecompress(batch, textOut);
        } else if (!gpuHostDecompress(batch, textOut)) {
            std::exit(1);
        }
        const auto t1 = Clock::now();
        if (measured) {
            nvtxRangeEnd(rng);
            cudaProfilerStop();
            times.push_back(std::chrono::duration<double>(t1 - t0).count());
        }
    }

    printf("{\"warmup\":%d,\"samples\":%d,\"batch\":%zu,\"batch_bytes\":%zu,\"times_s\":[",
           warmup, samples, preloaded.size(), batch_bytes);
    for (size_t j = 0; j < times.size(); ++j) {
        if (j) printf(",");
        printf("%.6f", times[j]);
    }
    printf("]}\n");
    fflush(stdout);
}

}  // namespace

int main(int argc, char** argv) {
    int iters = 5;
    std::string structure = "test/test.pdb";
    std::string sizesArg = "8,16,32,64,128,256";
    std::string inputDir;
    bool recursive = false;
    bool cpuMode = false;
    int warmup = 0;
    int benchSamples = 0;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&](const char* name) -> std::string {
            if (i + 1 >= argc) {
                fprintf(stderr, "[error] %s needs a value\n", name);
                std::exit(2);
            }
            return argv[++i];
        };
        if (a == "--iters") iters = parsePositiveInt("--iters", next("--iters"));
        else if (a == "--structure") structure = next("--structure");
        else if (a == "--sizes") sizesArg = next("--sizes");
        else if (a == "--max-batch-size")
            g_max_batch_size = parsePositiveInt("--max-batch-size", next("--max-batch-size"));
        else if (a == "--max-residues")
            g_max_residues = parsePositiveInt("--max-residues", next("--max-residues"));
        else if (a == "--input-dir") inputDir = next("--input-dir");
        else if (a == "--recursive") recursive = true;
        else if (a == "--cpu") cpuMode = true;
        else if (a == "--bench-warmup") warmup = parseNonNegativeInt("--bench-warmup", next("--bench-warmup"));
        else if (a == "--bench-samples") benchSamples = parsePositiveInt("--bench-samples", next("--bench-samples"));
        else if (a == "-h" || a == "--help") {
            printf("Usage: %s [--iters N] [--sizes 8,32,...] [--structure test/test.pdb]\n"
                   "       %s --input-dir DIR --bench-samples N [--bench-warmup N] [--recursive] [--cpu]\n",
                   argv[0], argv[0]);
            return 0;
        } else {
            fprintf(stderr, "[error] unknown arg: %s\n", a.c_str());
            return 2;
        }
    }

    std::string reason;
    if (!cudaRuntimeAvailable(&reason)) {
        fprintf(stderr, "CUDA runtime unavailable: %s\n", reason.c_str());
        return 1;
    }

    if (!inputDir.empty()) {
        if (benchSamples <= 0) {
            fprintf(stderr, "[error] --input-dir requires --bench-samples N\n");
            return 2;
        }
        runInputDirBench(inputDir, recursive, cpuMode, warmup, benchSamples);
        return 0;
    }

    const std::string base_fcz = toFcz(structure);
    const std::vector<int> sizes = parseSizes(sizesArg);

    printf("structure=%s  iters=%d  (median wall time)\n\n",
           structure.c_str(), iters);

    verifySingle(base_fcz);
    printf("\n");

    // Per-call GPU pipeline setup floor (empty batch) — amortized/warm, matches
    // the Python harness. Also warms the pipeline cache for the table below.
    std::vector<std::string> dummyOut;
    const double setup = timeit([&] { gpuHostDecompress({}, dummyOut); }, iters,
                                "setup floor (empty batch)");
    printf("GPU per-call pipeline setup floor (empty batch): %.2f ms "
           "- included in every GPU time below.\n\n", setup);

    printf("%5s | %8s | %9s | %10s | %9s | %7s | %9s\n",
           "N", "atoms", "CPU ms", "GPUhost ms", "GPUfmt ms",
           "CPU/GPU", "GPU str/s");
    printf("--------------------------------------------------------------------"
           "-------------------------\n");
    for (int n : sizes) {
        std::vector<std::string> batch(static_cast<size_t>(n), base_fcz);
        nvtx3::scoped_range r{"batch N=" + std::to_string(n)};
        const Row row = benchBatch(batch, iters);
        printRow(row);
    }

    printf("\nNotes:\n");
    printf("  * CPU     = decodeStructureToPDB per structure (CPU decompress -> PDB text).\n");
    printf("  * GPUhost = GPU pipeline + host writeSegmentsToPDB per structure.\n");
    printf("  * GPUfmt  = GPU pipeline with device PDB formatting + parallel line stitch.\n");
    printf("  * Same internal functions as the Python binding, minus the PyUnicode layer.\n");
    printf("  * GPU pipeline cache is warmed; the ~180ms first-call setup is amortized.\n");
    return 0;
}
