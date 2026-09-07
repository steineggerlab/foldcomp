#!/usr/bin/env python3
"""Benchmark CPU vs GPU Foldcomp decompression (no file output).

Measures, for a single structure and for batches of 8/16/32/64/128/256:

  * CPU decompress                 - foldcomp.decompress_batch(use_gpu=False)
  * GPU-host decompress            - foldcomp.decompress_batch(use_gpu=True), PDB text
  * GPU-formatted decompress       - foldcomp.decompress_batch(use_gpu=True, gpu_format=True)

Nothing is written to disk. Also verifies that single-file CPU and GPU
coordinates agree.

Run:  python test/benchmark_cpu_gpu.py [--iters N] [--structure test/test.pdb]
"""

import argparse
import ctypes
import statistics
import time
from pathlib import Path

import foldcomp

ROOT = Path(__file__).resolve().parent.parent

NVTX = None  # set in main()


def load_nvtx():
    """ctypes handle to libnvToolsExt for NVTX ranges, or None if unavailable.

    NVTX lets nsys attribute wall time to named ranges. It is optional: with no
    profiler attached the push/pop calls are cheap no-ops, and if the library is
    missing the benchmark still runs (ranges just don't appear in traces).
    """
    for name in ("libnvToolsExt.so.1", "libnvToolsExt.so"):
        try:
            lib = ctypes.CDLL(name)
        except OSError:
            continue
        lib.nvtxRangePushA.restype = ctypes.c_int
        lib.nvtxRangePushA.argtypes = [ctypes.c_char_p]
        lib.nvtxRangePop.restype = ctypes.c_int
        lib.nvtxRangePop.argtypes = []
        return lib
    return None


class nvtx_range:
    """Context manager pushing/popping a thread-local NVTX range.

    No-op when NVTX is unavailable, so it is safe to wrap any region.
    """

    __slots__ = ("_label",)

    def __init__(self, label):
        self._label = label.encode("ascii", "replace")

    def __enter__(self):
        if NVTX is not None:
            NVTX.nvtxRangePushA(self._label)
        return self

    def __exit__(self, *_exc):
        if NVTX is not None:
            NVTX.nvtxRangePop()
        return False


def pdb_coords_flat(results):
    """Flatten xyz from a list of (name, pdb_text) into a flat float list."""
    flat = []
    for _name, text in results:
        for line in text.splitlines():
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                flat.append(float(line[30:38]))
                flat.append(float(line[38:46]))
                flat.append(float(line[46:54]))
    return flat


def timeit(fn, iters, warmup=1, label=None):
    """Return the median wall time (seconds) of fn over `iters` runs.

    When `label` is given, each warmup/measured call is wrapped in an NVTX range
    (``<label>/warmup`` and ``<label>/iterK``) so a profiler shows exactly how
    much wall time each iteration of each case took.
    """
    tag = label or "timeit"
    for w in range(warmup):
        with nvtx_range(f"{tag}/warmup{w}"):
            fn()
    samples = []
    for k in range(iters):
        with nvtx_range(f"{tag}/iter{k}"):
            t0 = time.perf_counter()
            fn()
            samples.append(time.perf_counter() - t0)
    return statistics.median(samples)


def verify_single(base_fcz):
    """Verify single-file decompression.

    CPU vs GPU is only reported: Foldcomp's CPU and GPU reconstruction use
    different code paths and differ by ~0.2 A by design.
    """
    [(_, cpu_text)] = foldcomp.decompress_batch([base_fcz], use_gpu=False)
    [(_, gpu_text)] = foldcomp.decompress_batch([base_fcz], use_gpu=True)
    cpu = pdb_coords_flat([(None, cpu_text)])
    gpu_host = pdb_coords_flat([(None, gpu_text)])

    assert len(gpu_host) == len(cpu), "atom count mismatch"
    err_cpu_gpu = max(abs(cpu[i] - gpu_host[i]) for i in range(len(cpu)))
    print(f"Single-file: atoms={len(cpu) // 3}")
    print(
        f"  CPU     vs GPU      max|d|={err_cpu_gpu:.5f}  "
        f"(informational: CPU/GPU algorithms differ ~0.2 A)"
    )


def bench_batch(inputs, iters):
    """Return timing dict (seconds) for one batch of `inputs`."""
    n_struct = len(inputs)

    # --- CPU decompress (host PDB text; includes text formatting) ---
    def cpu_decompress():
        return foldcomp.decompress_batch(inputs, use_gpu=False)

    t_cpu = timeit(cpu_decompress, iters, label=f"CPU N={n_struct}")
    n_atoms = len(pdb_coords_flat(cpu_decompress())) // 3

    # --- GPU host decompress (use_gpu=True): GPU pipeline, but returns PDB text
    #     (device->host AtomCoordinate copy + text formatting on the host) ---
    def gpu_host_decompress():
        return foldcomp.decompress_batch(inputs, use_gpu=True)

    t_gpu_host = timeit(gpu_host_decompress, iters, label=f"GPUhost N={n_struct}")

    # --- GPU-formatted text (use_gpu=True, gpu_format=True): PDB ATOM records
    #     formatted on the device instead of single-threaded on the host ---
    def gpu_fmt_decompress():
        return foldcomp.decompress_batch(inputs, use_gpu=True, gpu_format=True)

    t_gpu_fmt = timeit(gpu_fmt_decompress, iters, label=f"GPUfmt N={n_struct}")

    return {
        "n_struct": n_struct,
        "n_atoms": n_atoms,
        "cpu": t_cpu,
        "gpu_host": t_gpu_host,
        "gpu_fmt": t_gpu_fmt,
    }


def fmt_row(r):
    def rate(t):
        return r["n_struct"] / t if t > 0 else float("inf")

    return (
        f"{r['n_struct']:>5} | {r['n_atoms']:>8} | "
        f"{r['cpu']*1e3:>9.2f} | {r['gpu_host']*1e3:>10.2f} | {r['gpu_fmt']*1e3:>9.2f} | "
        f"{r['cpu']/r['gpu_fmt']:>6.2f}x | {rate(r['gpu_fmt']):>8.0f}"
    )


def main():
    global NVTX
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--iters",
        type=int,
        default=5,
        help="measured iterations per case (median reported)",
    )
    parser.add_argument(
        "--structure",
        default=str(ROOT / "test" / "test.pdb"),
        help="PDB used as the repeated batch element",
    )
    parser.add_argument(
        "--sizes", default="8,16,32,64,128,256", help="comma-separated batch sizes"
    )
    args = parser.parse_args()

    if not foldcomp.cuda_available():
        raise SystemExit("CUDA runtime is unavailable; cannot benchmark GPU path")
    NVTX = load_nvtx()
    if NVTX is None:
        print(
            "[note] libnvToolsExt not found; NVTX ranges disabled "
            "(timings unaffected)\n"
        )

    pdb = Path(args.structure).read_bytes()
    base_fcz = foldcomp.compress(Path(args.structure).stem, pdb)
    sizes = [int(s) for s in args.sizes.split(",") if s.strip()]

    print(
        f"structure={args.structure}  iters={args.iters}  "
        f"(median wall time; GPU calls synchronize)\n"
    )

    with nvtx_range("verify_single"):
        verify_single(base_fcz)
    print()

    # Fixed per-call GPU pipeline setup floor (slot construction + buffer
    # preallocate + stream/event init), measured with an empty batch.
    setup = timeit(
        lambda: foldcomp.decompress_batch([], use_gpu=True),
        args.iters,
        label="setup floor (empty batch)",
    )
    print(
        f"GPU per-call pipeline setup floor (empty batch): {setup*1e3:.2f} ms "
        f"- included in every GPU time below.\n"
    )

    header = (
        f"{'N':>5} | {'atoms':>8} | "
        f"{'CPU ms':>9} | {'GPUhost ms':>10} | {'GPUfmt ms':>9} | "
        f"{'CPU/GPU':>7} | {'GPU str/s':>9}"
    )
    print(header)
    print("-" * len(header))
    for n in sizes:
        with nvtx_range(f"batch N={n}"):
            r = bench_batch([base_fcz] * n, args.iters)
        print(fmt_row(r))

    print("\nNotes:")
    print("  * CPU      = decompress_batch(use_gpu=False): CPU decompress -> PDB text.")
    print("  * GPUhost  = decompress_batch(use_gpu=True): GPU decompress, then")
    print("    device->host AtomCoordinate copy + PDB text formatting on the host.")
    print("  * GPUfmt   = decompress_batch(use_gpu=True, gpu_format=True): PDB ATOM")
    print("    records formatted on the GPU in-pipeline, no AtomCoordinate host copy")
    print("    (byte-identical output).")
    print("  * CPU/GPU speedup and str/s compare the CPU and GPUfmt paths.")
    print("  * The GPU pipeline (slots + pinned buffers) is cached across calls, so")
    print(
        "    the ~180ms first-call setup is amortized; foldcomp.gpu_pipeline_release()"
    )
    print("    frees it. Times above are warmed (cache hit).")


if __name__ == "__main__":
    main()
