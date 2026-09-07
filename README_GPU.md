# GPU Acceleration

Foldcomp has an optional CUDA backend that decompresses FCZ/FCMP/FCZC
structures on the GPU. It is off by default; the CPU path is unchanged and
remains the reference implementation. This document covers the GPU-specific
architecture, platform support, build prerequisites, the Python API, and a
known, intentional numerical divergence from the CPU reference.

## Architecture

The GPU path is a multi-threaded, batched pipeline (`src/gpu/gpu_decompression_pipeline.{cpp,h}`)
built around a pool of reusable pipeline instances:

- **Reader threads** parse incoming FCZ/FCZC entries in parallel and stage
  them into fixed-size batches.
- **Discretizer / continuize kernel** (`gpu_discretizer.cu`) reconstructs the
  quantized torsion angles on-device.
- **NeRF backbone kernel** (`gpu_nerf.cu`) reconstructs backbone atom
  coordinates per protein segment, forward and backward, blending at segment
  boundaries — the GPU counterpart of `Nerf::reconstructWithReversed` on the
  CPU.
- **Sidechain kernel** (`gpu_sidechain.cu`) places sidechain atoms from the
  reconstructed backbone.
- **Optional on-device PDB writer** (`gpu_pdb_writer.cu`) formats ATOM lines
  directly on the GPU, skipping a host-side text-assembly pass.

Each pipeline instance owns a fixed number of GPU "slots" (in-flight
sub-batches) plus pinned host staging buffers, and is checked out of a
process-wide pool (`pool_size`, default 2) keyed by config/mode/device — so
repeated calls (e.g. from a long-lived Python process) reuse allocations
instead of re-allocating every call. See `src/gpu/gpu_buffer.h` for the
underlying grow-only pinned/device buffer types.

Two output modes are supported at the pipeline level, selected by the
caller:

1. **Host segments** — atoms copied back to host memory as `AtomCoordinate`s
   (the default; used to write PDB/mmCIF). This is what the Python API below
   uses.
2. **GPU-formatted PDB lines** — ATOM lines formatted on-device, then copied
   back to a pinned host buffer as text, avoiding host-side formatting
   entirely (`decompress_batch(..., gpu_format=True)`).

A third mode, device-resident coordinates (decompressed atoms left as a
`float32 [N, 3]` buffer in GPU memory, never copied back to the host), existed
in the pipeline but was never wired up to the Python API — it backed the
since-removed `foldcomp.decompress_batch_gpu()` / `__cuda_array_interface__`
surface — and has been removed as unreachable code.

## Supported platforms

- **GPU decompression**: Linux and Windows, with an NVIDIA GPU and a working
  CUDA driver/toolkit. macOS is not supported for the GPU path (no CUDA on
  macOS).
- **CPU decompression**: unchanged and available everywhere Foldcomp already
  builds (Linux, Windows/MSVC, macOS).
- If no usable CUDA runtime/device is detected, Foldcomp automatically falls
  back to the CPU path — `foldcomp.decompress_batch(..., use_gpu=None)`
  (the default) selects GPU only when available, and the CLI's `--gpu` flag
  (the default for CUDA-enabled builds) falls back to `--cpu` with a warning.

## Build prerequisites

- CMake **3.18** or newer.
- A CUDA Toolkit with an `nvcc` (or Clang) that supports C++17 device code
  (CUDA **12.0+** recommended) — required for `find_package(CUDAToolkit)` and
  `enable_language(CUDA)` in `cmake/CudaBuild.cmake`.
- An NVIDIA GPU of compute capability 6.0 (Pascal) or newer. Target
  architectures are auto-discovered by probing the active compiler against a
  candidate list (`sm_60` through `sm_120`); pass
  `-DCMAKE_CUDA_ARCHITECTURES=...` explicitly to override discovery. If no
  candidate is accepted, the build falls back to `sm_80` with a warning.
- OpenMP (already required for the CPU build's library/Python targets).

Enable the GPU build with:

```bash
cmake -DBUILD_CUDA=ON ..
```

or, for a local Python build via scikit-build:

```bash
CMAKE_ARGS="-DBUILD_CUDA=ON" pip install .
```

A CPU-only build (`-DBUILD_CUDA=OFF`, the default) compiles and links with no
CUDA dependency at all.

## Executable

GPU decompression is used automatically by `foldcomp decompress` in
CUDA-enabled builds; the options below (a subset of `foldcomp --help`, see
the main [README](README.md#executable) for the rest) control it:

```
 --cpu                    use original single-threaded CPU decompression path (no GPU)
 --gpu                    use GPU decompression path [default for CUDA-enabled builds;
                          falls back to --cpu with a warning if no CUDA runtime is available]
 --memory-only            decompress to memory only, skip disk/tar/db writes (GPU decompress path only)
 --max-batch-size N       max structures per GPU batch [default=100]
 --max-residues N         upper bound on residues per structure for pinned buffer pre-allocate [default=3000]
 --reader-threads N       parallel FCZ parser threads for GPU batch decompression [default=<nproc>]
 --no-fused-disc          disable fused continuize kernels for discretizer (GPU only) [default=false]
 --sort-by-length         pre-sort input files by protein length (nResidue) for uniform GPU batches
```

## Python API

CUDA-enabled Python builds expose GPU decompression through the same
`decompress_batch` entry point as the CPU path, plus a couple of GPU-specific
helpers:

```python
import foldcomp

with open("example.fcz", "rb") as fh:
    fcz_bytes = fh.read()

if foldcomp.cuda_available():
    pdb_texts = foldcomp.decompress_batch(
        [fcz_bytes],       # list of FCZ/FCMP/FCZC byte strings
        use_gpu=True,      # None (default): GPU if available, else CPU
        format="pdb",      # or "mmcif"
        gpu_format=True,   # format ATOM lines on-device (mode 3 above)
        pool_size=2,       # concurrent pipeline instances kept warm
    )
```

**Resource allocation.** The first `decompress_batch(use_gpu=True, ...)` call
for a given `(max_batch_size, max_residues, gpu_slots, pool_size, mode,
device)` configuration builds (or checks out) a pooled pipeline instance —
GPU device buffers, pinned host staging buffers, and CUDA
streams/events — sized up front. That instance is then **cached and reused**
across subsequent calls instead of being freed, so repeated calls from a
long-lived process (e.g. a service handling many requests) don't pay
allocation cost every time. The tradeoff is that this memory stays resident
until the process exits or is explicitly released:

```python
foldcomp.gpu_pipeline_release()   # free all cached pipeline instances now

# or, scope it to a block:
with foldcomp.gpu_context():
    pdb_texts = foldcomp.decompress_batch(fcz_list, use_gpu=True)
# gpu_pipeline_release() runs automatically here; the next
# decompress_batch(use_gpu=True) call rebuilds the pool on demand.
```

`gpu_context()` is a single global resource, not scoped per-call — don't nest
it or use it concurrently from multiple threads, since exiting one block
releases pipelines that another still-open block or thread may depend on.

## Known accuracy divergence from the CPU reference

GPU decompression is designed to closely match the CPU reference output, but
it is **not bit-exact**. Each segment's backbone is reconstructed in three
passes:

1. **Forward** — walk the segment from its start anchor.
2. **Backward + blend** — walk it from its end anchor and blend with the
   forward pass (weight `m/total` at local atom `m`).
3. **Rigid-frame-transfer correction** (`backbone_segment_frame_transfer_kernel`,
   `src/gpu/gpu_nerf.cu`) — re-expresses each segment's interior atoms in the
   coordinate frame of the *previous* segment's actual blended output, instead
   of the exact stored anchor the forward/backward passes used to seed from.

On the CPU, each segment is seeded from the CPU's own *blended* result for
the previous segment's boundary residue (`src/foldcomp.cpp:1049-1052`) — a
serial, segment-to-segment dependency. The GPU's forward/backward passes seed
from the exact stored anchor instead, so segments can be reconstructed
independently and in parallel; pass 3 then corrects that seed choice after
the fact by rigidly transforming the result into the frame the CPU's blended
seed would have produced, without reintroducing the serial dependency. (Pass
3 also does not need to correct the segment's own start/end boundary
residues, which are otherwise already exact or handled by the blend — see
the kernel's doc comment for the exact indexing.)

**Measured residual** (`test/test.pdb`, GPU vs. CPU, same binary — reproduce
with `test/compare_rmsd_dirs.py`, decompressing the same input with `--cpu`
and `--gpu` into two directories and diffing them):

| | max per-atom deviation (GPU vs. CPU) | mean |
|---|---|---|
| as shipped (pass 3 enabled) | 0.0064 Å | 0.0008 Å |
| pass 3 disabled | 0.0203 Å | 0.0057 Å |

**Pass 3 improves agreement with the CPU reference — it does not improve
accuracy against the original input structure**, and should not be "optimised
away" on the assumption that it does. Backbone RMSD against the original
`test/test.pdb` barely moves either way: CPU reference 0.045377 Å, GPU with
pass 3 (as shipped) 0.045548 Å, GPU with pass 3 disabled 0.045672 Å. Pass 3
exists purely to make the GPU path a closer drop-in match for the CPU path's
output, not to reconstruct the original structure more faithfully.

This residual is small relative to the format's overall reconstruction error
and is exercised by the tolerance bounds in `test/test_e2e_gpu.sh`.
