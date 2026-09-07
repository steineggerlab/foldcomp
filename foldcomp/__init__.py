import contextlib

from . import foldcomp as _foldcomp
from .foldcomp import *
from .setup import setup, setup_async
from .util import split_pdb_by_chain


@contextlib.contextmanager
def gpu_context():
    """Release cached GPU decompression resources when the block exits.

    decompress_batch(use_gpu=True) caches pooled GPU pipeline instances
    (device/pinned buffers, CUDA streams) for reuse across
    calls, freed only on process exit or an explicit gpu_pipeline_release().
    Wrapping a block in `with foldcomp.gpu_context():` calls
    gpu_pipeline_release() automatically once the block exits, so GPU memory
    isn't held past the point where it's needed.

    The underlying pool is a single global resource (not scoped to this
    context), so this is safe for sequential, single-threaded use only. Do
    not nest gpu_context() blocks or use them concurrently from multiple
    threads: exiting one releases pipelines that another still-open block or
    thread may be relying on staying warm.
    """
    try:
        yield
    finally:
        if hasattr(_foldcomp, "gpu_pipeline_release"):
            _foldcomp.gpu_pipeline_release()
