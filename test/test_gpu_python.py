from pathlib import Path

import foldcomp
import pytest

ROOT = Path(__file__).resolve().parent.parent


def _multichain_fcz():
    pdb = (ROOT / "test" / "multichain.pdb").read_bytes()
    return foldcomp.compress("multi", pdb)


def _require_gpu():
    if not foldcomp.cuda_available():
        pytest.skip("CUDA runtime is unavailable")


def test_cuda_available_is_runtime_boolean():
    assert isinstance(foldcomp.cuda_available(), bool)


def test_decompress_batch_cpu_preserves_order_and_chains():
    multi = _multichain_fcz()
    single = foldcomp.compress("single", (ROOT / "test" / "test.pdb").read_bytes())
    results = foldcomp.decompress_batch([multi, single], use_gpu=False)
    assert [name for name, _ in results] == ["multi", "single"]
    chains = {
        line[21]
        for line in results[0][1].splitlines()
        if line.startswith(("ATOM", "HETATM")) and len(line) > 21
    }
    assert chains == {"A", "B"}


def test_decompress_batch_cpu_supports_mmcif():
    [(name, text)] = foldcomp.decompress_batch(
        [_multichain_fcz()], use_gpu=False, format="mmcif"
    )
    assert name == "multi"
    assert text.startswith("data_multi")
    assert "_atom_site.auth_asym_id" in text


def test_decompress_batch_auto_mode_always_produces_output():
    [(name, text)] = foldcomp.decompress_batch([_multichain_fcz()], use_gpu=None)
    assert name == "multi"
    assert any(line.startswith("ATOM") for line in text.splitlines())


def test_decompress_batch_gpu_handles_container_fragments_across_batches():
    _require_gpu()
    multichain = _multichain_fcz()
    [(name, text)] = foldcomp.decompress_batch(
        [multichain], use_gpu=True, max_batch_size=1
    )
    assert name == "multi"
    chains = {
        line[21]
        for line in text.splitlines()
        if line.startswith(("ATOM", "HETATM")) and len(line) > 21
    }
    assert chains == {"A", "B"}
    [(_, cpu_text)] = foldcomp.decompress_batch([multichain], use_gpu=False)

    def identities(pdb_text):
        return [
            (line[12:16], line[17:20], line[21:22], line[22:26], line[26:27])
            for line in pdb_text.splitlines()
            if line.startswith(("ATOM", "HETATM"))
        ]

    assert identities(text) == identities(cpu_text)


def test_decompress_batch_gpu_handles_raw_container_fragments():
    _require_gpu()
    raw_container = foldcomp.compress(
        "raw", (ROOT / "test" / "test.pdb").read_bytes(), max_backbone_rmsd=0.0
    )
    [(name, text)] = foldcomp.decompress_batch([raw_container], use_gpu=True)
    assert name == "raw"
    assert any(line.startswith("ATOM") for line in text.splitlines())


def test_decompress_batch_gpu_preserves_model_and_oxt():
    _require_gpu()
    source_lines = (ROOT / "test" / "test_af.pdb").read_text().splitlines()
    shifted_atoms = []
    for line in source_lines:
        if not line.startswith(("ATOM", "HETATM", "TER")):
            continue
        residue_index = int(line[22:26]) + 41
        shifted_atoms.append(f"{line[:22]}{residue_index:4d}{line[26:]}")
    atoms = "\n".join(shifted_atoms)
    two_models = (
        f"MODEL        1\n{atoms}\nENDMDL\nMODEL        2\n{atoms}\nENDMDL\nEND\n"
    )
    model_container = foldcomp.compress("model_oxt", two_models.encode())
    [(name, text)] = foldcomp.decompress_batch([model_container], use_gpu=True)
    assert name == "model_oxt"
    assert "MODEL" in text
    assert any(" OXT " in line for line in text.splitlines())
    [(_, cpu_text)] = foldcomp.decompress_batch([model_container], use_gpu=False)
    gpu_identities = [
        line[12:27] for line in text.splitlines() if line.startswith(("ATOM", "HETATM"))
    ]
    cpu_identities = [
        line[12:27]
        for line in cpu_text.splitlines()
        if line.startswith(("ATOM", "HETATM"))
    ]
    assert gpu_identities == cpu_identities


def test_decompress_batch_rejects_invalid_inputs():
    with pytest.raises(TypeError):
        foldcomp.decompress_batch(b"not-a-list")
    with pytest.raises(TypeError):
        foldcomp.decompress_batch(["not-bytes"])
    with pytest.raises(TypeError):
        foldcomp.decompress_batch([], use_gpu=1)
    with pytest.raises(ValueError, match="format"):
        foldcomp.decompress_batch([], format="xyz")


def test_explicit_gpu_fails_when_runtime_is_unavailable():
    if foldcomp.cuda_available():
        pytest.skip("CUDA runtime is available")
    with pytest.raises(RuntimeError, match="CUDA"):
        foldcomp.decompress_batch([_multichain_fcz()], use_gpu=True)


def _pdb_xyz(text):
    return [
        (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        for line in text.splitlines()
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 54
    ]


def _raw_container():
    raw = foldcomp.compress(
        "raw", (ROOT / "test" / "test.pdb").read_bytes(), max_backbone_rmsd=0.0
    )
    assert raw.startswith(b"FCZC")
    return raw


# Two adjacent chains ("AA", "AB") that share their first character. PDB format
# only has one column for chain, so both render as "A" there -- but the segment
# boundary between them (TER placement) must still be decided from the full
# auth_asym_id, not just that shared first byte. Requires mmCIF input since PDB
# format itself can't carry a >1-char chain ID.
_WIDE_CHAIN_CIF = b"""\
data_TWO_CHAIN_TEST
#
_entry.id 'TWO_CHAIN_TEST'
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . GLY AA 1 1 ? 0.000 0.000 0.000 1 20.0 ? 1 AA 1
ATOM 2 CA CA . GLY AA 1 1 ? 1.458 0.000 0.000 1 20.0 ? 1 AA 1
ATOM 3 C C . GLY AA 1 1 ? 2.009 1.420 0.000 1 20.0 ? 1 AA 1
ATOM 4 O O . GLY AA 1 1 ? 1.251 2.390 0.000 1 20.0 ? 1 AA 1
ATOM 5 N N . GLY AB 1 1 ? 10.000 10.000 0.000 1 20.0 ? 1 AB 1
ATOM 6 CA CA . GLY AB 1 1 ? 11.458 10.000 0.000 1 20.0 ? 1 AB 1
ATOM 7 C C . GLY AB 1 1 ? 12.009 11.420 0.000 1 20.0 ? 1 AB 1
ATOM 8 O O . GLY AB 1 1 ? 11.251 12.390 0.000 1 20.0 ? 1 AB 1
#
"""


def _wide_chain_container():
    # max_backbone_rmsd=0.0 forces raw-atom container fragments (like
    # _raw_container above), which is where the wide ChainId is carried.
    container = foldcomp.compress(
        "wide_chain", _WIDE_CHAIN_CIF, format="cif", max_backbone_rmsd=0.0
    )
    assert container.startswith(b"FCZC")
    return container


# One coordinate (12345.678, the reviewer's own example) that overflows the
# 8-char fixed-width PDB coordinate field. max_backbone_rmsd=0.0 forces a
# raw-atom container fragment, so the input float passes through both
# backends unmodified (no NeRF reconstruction to round it away).
_OVERFLOW_COORD_CIF = b"""\
data_OVERFLOW_TEST
#
_entry.id 'OVERFLOW_TEST'
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . GLY A 1 1 ? 12345.678 0.000 0.000 1 20.0 ? 1 A 1
ATOM 2 CA CA . GLY A 1 1 ? 1.458 0.000 0.000 1 20.0 ? 1 A 1
ATOM 3 C C . GLY A 1 1 ? 2.009 1.420 0.000 1 20.0 ? 1 A 1
ATOM 4 O O . GLY A 1 1 ? 1.251 2.390 0.000 1 20.0 ? 1 A 1
#
"""


def _overflow_coord_container():
    container = foldcomp.compress(
        "overflow", _OVERFLOW_COORD_CIF, format="cif", max_backbone_rmsd=0.0
    )
    assert container.startswith(b"FCZC")
    return container


def test_pdb_coordinate_overflow_clamped_on_cpu():
    # Backward compatibility / regression: appendRightAlignedFixed (host)
    # used to widen an out-of-range coordinate field instead of clamping it
    # like the GPU formatter's pf_fixed does (e.g. "12345.679" instead of
    # "********"), silently shifting every following column. It must now
    # clamp on the CPU path too, matching the device formatter.
    container = _overflow_coord_container()
    [(_, pdb)] = foldcomp.decompress_batch([container], use_gpu=False)
    atom_line = next(line for line in pdb.splitlines() if line.startswith("ATOM"))
    assert "********" in atom_line
    assert "12345" not in atom_line


def test_pdb_coordinate_overflow_matches_between_cpu_and_gpu():
    _require_gpu()
    container = _overflow_coord_container()
    [(_, cpu_pdb)] = foldcomp.decompress_batch([container], use_gpu=False)
    [(_, gpu_host_pdb)] = foldcomp.decompress_batch([container], use_gpu=True)
    [(_, gpu_line_pdb)] = foldcomp.decompress_batch(
        [container], use_gpu=True, gpu_format=True
    )
    assert gpu_host_pdb == cpu_pdb
    assert gpu_line_pdb == cpu_pdb


def test_gpu_pdb_writer_byte_exact():
    _require_gpu()
    single = foldcomp.compress("single", (ROOT / "test" / "test.pdb").read_bytes())
    raw = _raw_container()
    inputs = [single, _multichain_fcz(), raw]
    gpu_host_formatted = foldcomp.decompress_batch(inputs, use_gpu=True)
    gpu = foldcomp.decompress_batch(inputs, use_gpu=True, gpu_format=True)
    assert [n for n, _ in gpu_host_formatted] == [n for n, _ in gpu]
    for (_, ht), (_, gt) in zip(gpu_host_formatted, gpu):
        assert ht == gt  # GPU-formatted PDB is byte-identical to the host formatter


def test_wide_chain_id_ter_boundary_cpu():
    # Backward compatibility: the CPU path (writeSegmentsToPDB) already compared
    # the full chain id and must keep separating "AA"/"AB" with TER.
    container = _wide_chain_container()
    [(_, pdb)] = foldcomp.decompress_batch([container], use_gpu=False)
    assert pdb.count("TER") == 2
    lines = pdb.splitlines()
    ter_lines = [line for line in lines if line.startswith("TER")]
    assert len(ter_lines) == 2


def test_wide_chain_id_ter_boundary_gpu_host_formatted():
    _require_gpu()
    container = _wide_chain_container()
    [(_, pdb)] = foldcomp.decompress_batch([container], use_gpu=True)
    assert pdb.count("TER") == 2


def test_wide_chain_id_ter_boundary_gpu_line_formatted():
    # Regression for the GPU-formatted (gpu_format=True) line-driven PDB
    # assembler: PdbLineSegment::chain used to be a single char, so segments
    # whose full chain ids shared only their first character (here "AA"/"AB",
    # which both render as "A" in the single PDB chain column) were treated as
    # the same chain and the TER between them was silently dropped.
    _require_gpu()
    container = _wide_chain_container()
    [(_, pdb)] = foldcomp.decompress_batch([container], use_gpu=True, gpu_format=True)
    assert pdb.count("TER") == 2


def test_wide_chain_id_gpu_line_formatted_matches_cpu_and_host_formatted():
    _require_gpu()
    container = _wide_chain_container()
    [(_, cpu_pdb)] = foldcomp.decompress_batch([container], use_gpu=False)
    [(_, gpu_host_pdb)] = foldcomp.decompress_batch([container], use_gpu=True)
    [(_, gpu_line_pdb)] = foldcomp.decompress_batch(
        [container], use_gpu=True, gpu_format=True
    )
    assert gpu_host_pdb == cpu_pdb
    assert gpu_line_pdb == cpu_pdb
    assert gpu_line_pdb == gpu_host_pdb


def test_gpu_pdb_writer_ignored_for_non_pdb_and_cpu():
    _require_gpu()
    fcz = foldcomp.compress("t", (ROOT / "test" / "test.pdb").read_bytes())
    # mmcif ignores gpu_format (still host mmcif)
    [(_, cif)] = foldcomp.decompress_batch(
        [fcz], use_gpu=True, format="mmcif", gpu_format=True
    )
    assert cif.startswith("data_")
    # CPU path ignores gpu_format
    [(_, a)] = foldcomp.decompress_batch([fcz], use_gpu=False, gpu_format=True)
    [(_, b)] = foldcomp.decompress_batch([fcz], use_gpu=False)
    assert a == b


def test_decompress_batch_gpu_pipeline_cache_and_release():
    _require_gpu()
    multi = _multichain_fcz()

    def coords():
        return _pdb_xyz(foldcomp.decompress_batch([multi], use_gpu=True)[0][1])

    c1 = coords()  # builds + caches the pipeline
    c2 = coords()  # reuses the cache
    foldcomp.gpu_pipeline_release()  # frees the cache
    c3 = coords()  # rebuilds the cache

    assert c1 == c2 == c3  # identical results across cache/rebuild
    assert foldcomp.gpu_pipeline_release() is None  # idempotent, returns None


def test_gpu_context_releases_pipeline_on_exit():
    _require_gpu()
    multi = _multichain_fcz()

    def coords():
        return _pdb_xyz(foldcomp.decompress_batch([multi], use_gpu=True)[0][1])

    with foldcomp.gpu_context():
        c1 = coords()  # builds + caches the pipeline
        c2 = coords()  # reuses the cache, still inside the block
    # Block exited: gpu_context() should have released the cache, same as an
    # explicit gpu_pipeline_release() call — the next use just rebuilds it.
    c3 = coords()

    assert c1 == c2 == c3


def test_gpu_context_releases_on_exception():
    _require_gpu()
    with pytest.raises(ValueError):
        with foldcomp.gpu_context():
            foldcomp.decompress_batch([_multichain_fcz()], use_gpu=True)
            raise ValueError("boom")
    # Cache was still released despite the exception; a fresh call works.
    [(_, text)] = foldcomp.decompress_batch([_multichain_fcz()], use_gpu=True)
    assert len(_pdb_xyz(text)) > 0


# --------------------------------------------------------------------------
# Pipeline pool: concurrent decompress_batch(use_gpu=True) calls from multiple
# Python threads (GIL is released around the C++ call — see foldcomp.cxx).
# These exercise the pooled-pipeline concurrency (each concurrent caller gets
# its own GPU slots + batch pool + capture buffer instead of fully
# serializing on one lock).
# --------------------------------------------------------------------------


def test_decompress_batch_gpu_concurrent_calls_are_correct():
    _require_gpu()
    import concurrent.futures

    inputs = [
        foldcomp.compress("single", (ROOT / "test" / "test.pdb").read_bytes()),
        _multichain_fcz(),
        _raw_container(),
    ]
    # Reference: single-threaded PDB text per input.
    refs = [
        _pdb_xyz(text) for _, text in foldcomp.decompress_batch(inputs, use_gpu=True)
    ]

    def run_one(inp, ref):
        [(_, text)] = foldcomp.decompress_batch([inp], use_gpu=True, pool_size=4)
        dev = _pdb_xyz(text)
        assert len(dev) == len(ref)
        for k in range(len(dev)):
            for c in range(3):
                assert abs(dev[k][c] - ref[k][c]) < 0.01

    # Each thread repeatedly decompresses a distinct input, so cross-thread
    # corruption (e.g. two runs sharing a capture buffer) would show up as a
    # coordinate mismatch against that input's own single-threaded reference.
    with concurrent.futures.ThreadPoolExecutor(max_workers=len(inputs)) as pool:
        futures = [
            pool.submit(run_one, inp, ref)
            for _ in range(5)
            for inp, ref in zip(inputs, refs)
        ]
        for f in futures:
            f.result()


def test_decompress_batch_gpu_concurrent_calls_dont_serialize():
    _require_gpu()
    import concurrent.futures
    import threading
    import time

    # A raw wall-clock throughput comparison (N concurrent vs N serial calls) is
    # too hardware-dependent to be a reliable test: on a single, already
    # near-saturated GPU, oversubscribing it with concurrent kernels can be
    # *slower* in aggregate than running them back-to-back, even though the
    # pool is working correctly. Instead, test the property the pool actually
    # guarantees: a small concurrent call is not blocked behind another
    # thread's much larger in-flight run. Under the old single-global-cache
    # design, the small call couldn't even start its run until the big call's
    # entire run finished (one mutex held for the whole call), so it would
    # always finish after. With pool_size=2, the small call gets its own
    # pooled instance and finishes independently.
    multi = _multichain_fcz()
    big_batch = [multi] * 8000
    small_batch = [multi]

    # Warm up both configs' pooled instances outside the timed region.
    foldcomp.decompress_batch(big_batch, use_gpu=True, pool_size=2)
    foldcomp.decompress_batch(small_batch, use_gpu=True, pool_size=2)

    big_started = threading.Event()
    timings = {}

    def run_big():
        big_started.set()
        foldcomp.decompress_batch(big_batch, use_gpu=True, pool_size=2)
        timings["big_end"] = time.perf_counter()

    def run_small():
        big_started.wait()
        time.sleep(0.05)  # let the big call's run actually get underway
        foldcomp.decompress_batch(small_batch, use_gpu=True, pool_size=2)
        timings["small_end"] = time.perf_counter()

    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        f_big = pool.submit(run_big)
        f_small = pool.submit(run_small)
        f_big.result()
        f_small.result()

    assert timings["small_end"] < timings["big_end"], (
        "small concurrent call finished after the big call — "
        "looks like calls are still fully serializing"
    )


def test_decompress_batch_gpu_concurrent_calls_pool_size_one():
    _require_gpu()
    import concurrent.futures

    multi = _multichain_fcz()
    ref = _pdb_xyz(foldcomp.decompress_batch([multi], use_gpu=True)[0][1])

    def run_one():
        [(_, text)] = foldcomp.decompress_batch([multi], use_gpu=True, pool_size=1)
        dev = _pdb_xyz(text)
        assert len(dev) == len(ref)
        for k in range(len(dev)):
            for c in range(3):
                assert abs(dev[k][c] - ref[k][c]) < 0.01

    # pool_size=1 forces every concurrent call onto the same instance; correctness
    # must hold via the blocking checkout path, not just via having spare capacity.
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(run_one) for _ in range(8)]
        for f in futures:
            f.result()


def test_decompress_batch_gpu_pool_eviction_across_configs():
    _require_gpu()
    multi = _multichain_fcz()
    ref = _pdb_xyz(foldcomp.decompress_batch([multi], use_gpu=True)[0][1])

    # pool_size=1 forces every call to a different max_batch_size to evict and
    # rebuild the sole resident instance; correctness must survive that churn.
    for max_batch_size in (100, 1, 100, 1):
        [(_, text)] = foldcomp.decompress_batch(
            [multi], use_gpu=True, max_batch_size=max_batch_size, pool_size=1
        )
        dev = _pdb_xyz(text)
        assert len(dev) == len(ref)
        for k in range(len(dev)):
            for c in range(3):
                assert abs(dev[k][c] - ref[k][c]) < 0.01
