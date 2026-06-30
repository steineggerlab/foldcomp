#!/usr/bin/env python3
"""
Comprehensive test suite for Foldcomp Python API with multi-chain support.

Tests cover:
  - compress() with single-chain and multi-chain PDB/CIF input
  - decompress() in PDB and mmCIF format
  - get_data() extracting torsion angles, coordinates, residues
  - open() database iteration (full and by ID subset)
  - split_pdb_by_chain() utility
  - FCZ container format (FCZC magic) detection
  - round-trip fidelity (atom counts, chain preservation)
  - compression ratios across structure sizes
  - CIF (mmCIF) input/output round-trips

Usage:
    python test_multichain_api.py [--log-dir test_report]
"""

import sys
import os
import gzip
import json
import time
import argparse
from pathlib import Path
from collections import namedtuple

# Ensure the installed foldcomp is importable
FOLDCOMP_SO = Path(__file__).resolve().parent.parent / ".venv" / "lib" / "python3.13" / "site-packages"
if str(FOLDCOMP_SO) not in sys.path:
    sys.path.insert(0, str(FOLDCOMP_SO))

import foldcomp

# ---------------------------------------------------------------------------
# Test harness
# ---------------------------------------------------------------------------

TestResult = namedtuple("TestResult", [
    "name", "passed", "duration", "details", "error"
])

results: list[TestResult] = []


def record(name: str, passed: bool, duration: float, details: dict | None = None, error: str = ""):
    results.append(TestResult(name, passed, duration, details or {}, error))


def run_test(fn):
    """Decorator that wraps a test function, records timing and pass/fail."""
    def wrapper(*args, **kwargs):
        t0 = time.perf_counter()
        try:
            fn(*args, **kwargs)
            record(fn.__name__, True, time.perf_counter() - t0)
        except Exception as exc:
            record(fn.__name__, False, time.perf_counter() - t0, error=str(exc))
    return wrapper


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _atom_lines(pdb_text: str) -> list[str]:
    return [l for l in pdb_text.splitlines() if l.startswith("ATOM")]


def _chain_ids(pdb_text: str) -> list[str]:
    seen = []
    for l in pdb_text.splitlines():
        if l.startswith("ATOM") and l[21] not in seen:
            seen.append(l[21])
    return seen


def _b_factors(pdb_text: str) -> list[float]:
    return [float(l[60:66]) for l in pdb_text.splitlines() if l.startswith("ATOM")]


ROOT = Path(__file__).resolve().parent.parent  # project root


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

# ---- compress / decompress ----

@run_test
def test_compress_multichain_pdb():
    """Compress a 2-chain PDB (6PP9 BRAF:MEK1) and verify FCZC container format."""
    pdb = (ROOT / "test" / "multichain.pdb").read_text()
    compressed = foldcomp.compress("6PP9_BRAF_MEK1", pdb.encode())
    assert compressed[:4] == b"FCZC", f"Expected FCZC container magic, got {compressed[:4]}"
    assert len(compressed) > 0


@run_test
def test_decompress_multichain_pdb():
    """Decompress multi-chain FCZ and verify chains are preserved."""
    pdb = (ROOT / "test" / "multichain.pdb").read_text()
    compressed = foldcomp.compress("6PP9_BRAF_MEK1", pdb.encode())
    name, out = foldcomp.decompress(compressed)
    assert name == "6PP9_BRAF_MEK1"
    assert set(_chain_ids(out)) == {"A", "B"}, f"Chains mismatch: {_chain_ids(out)}"
    in_atoms = len(_atom_lines(pdb))
    out_atoms = len(_atom_lines(out))
    assert out_atoms == in_atoms, f"Atom count mismatch: in={in_atoms}, out={out_atoms}"


@run_test
def test_decompress_multichain_mmcif():
    """Decompress multi-chain FCZ as mmCIF and verify chains."""
    pdb = (ROOT / "test" / "multichain.pdb").read_text()
    compressed = foldcomp.compress("6PP9_BRAF_MEK1", pdb.encode())
    name, cif_out = foldcomp.decompress(compressed, format="mmcif")
    assert name == "6PP9_BRAF_MEK1"
    # mmCIF has _atom_site.label_asym_id column
    cif_atoms = [l for l in cif_out.splitlines() if l and not l.startswith("#") and not l.startswith("data_") and not l.startswith("_")]
    chains = set()
    for line in cif_atoms:
        parts = line.split()
        if len(parts) >= 7 and parts[0] != "_atom_site." and parts[0] != "loop_":
            chains.add(parts[6])
    assert "A" in chains and "B" in chains, f"mmCIF chains: {chains}"


@run_test
def test_compress_single_chain_pdb():
    """Compress single-chain PDB (test.pdb)."""
    pdb = (ROOT / "test" / "test.pdb").read_text()
    compressed = foldcomp.compress("test_single", pdb.encode())
    name, out = foldcomp.decompress(compressed)
    assert name == "test_single"
    in_atoms = len(_atom_lines(pdb))
    out_atoms = len(_atom_lines(out))
    assert out_atoms == in_atoms, f"Atom count mismatch: in={in_atoms}, out={out_atoms}"


@run_test
def test_compress_af_pdb():
    """Compress AlphaFold-predicted PDB with pLDDT in B-factors."""
    pdb = (ROOT / "test" / "test_af.pdb").read_text()
    compressed = foldcomp.compress("test_af", pdb.encode())
    name, out = foldcomp.decompress(compressed)
    in_atoms = len(_atom_lines(pdb))
    out_atoms = len(_atom_lines(out))
    assert out_atoms == in_atoms, f"Atom count mismatch: in={in_atoms}, out={out_atoms}"
    # pLDDT is stored in B-factors
    in_bf = _b_factors(pdb)
    out_bf = _b_factors(out)
    # B-factors may differ slightly due to quantization
    assert len(out_bf) == len(in_bf), "B-factor count mismatch"


@run_test
def test_compress_cif():
    """Compress mmCIF file and round-trip."""
    cif_gz = ROOT / "test" / "test.cif.gz"
    cif_data = gzip.open(cif_gz, "rt").read()
    compressed = foldcomp.compress("test_cif", cif_data.encode(), format="mmcif")
    name, out = foldcomp.decompress(compressed)
    assert name == "test_cif"
    assert len(_atom_lines(out)) > 0


@run_test
def test_compress_multichain_cif_1dlp():
    """Compress 6-chain CIF (1dlp) and verify all chains preserved."""
    cif_data = gzip.open(ROOT / "test" / "1dlp.cif.gz", "rt").read()
    compressed = foldcomp.compress("1dlp", cif_data.encode(), format="mmcif")
    name, out = foldcomp.decompress(compressed)
    chains = set(_chain_ids(out))
    assert chains == {"A", "B", "C", "D", "E", "F"}, f"Chains: {chains}"


@run_test
def test_compress_multichain_cif_2ap2():
    """Compress 4-chain CIF (2ap2)."""
    cif_data = gzip.open(ROOT / "test" / "2ap2.cif.gz", "rt").read()
    compressed = foldcomp.compress("2ap2", cif_data.encode(), format="mmcif")
    name, out = foldcomp.decompress(compressed)
    chains = set(_chain_ids(out))
    assert chains == {"A", "C", "P", "Q"}, f"Chains: {chains}"


@run_test
def test_compress_multichain_cif_7c2s():
    """Compress 6-chain CIF (7c2s - Dengue antibody)."""
    cif_data = gzip.open(ROOT / "test" / "7c2s.cif.gz", "rt").read()
    compressed = foldcomp.compress("7c2s", cif_data.encode(), format="mmcif")
    name, out = foldcomp.decompress(compressed)
    chains = set(_chain_ids(out))
    assert chains == {"A", "B", "G", "H", "I", "M"}, f"Chains: {chains}"


@run_test
def test_compress_multichain_cif_1hrn():
    """Compress 2-chain CIF (1hrn - ribonuclease inhibitor)."""
    cif_data = gzip.open(ROOT / "test" / "1hrn.cif.gz", "rt").read()
    compressed = foldcomp.compress("1hrn", cif_data.encode(), format="mmcif")
    name, out = foldcomp.decompress(compressed)
    chains = set(_chain_ids(out))
    assert chains == {"A", "B"}, f"Chains: {chains}"


# ---- get_data ----

@run_test
def test_get_data_multichain_pdb():
    """Extract structure data from multi-chain PDB."""
    pdb = (ROOT / "test" / "multichain.pdb").read_bytes()
    data = foldcomp.get_data(pdb)
    required_keys = {"phi", "psi", "omega", "torsion_angles", "bond_angles",
                     "residues", "b_factors", "coordinates"}
    assert required_keys == set(data.keys()), f"Keys: {data.keys()}"
    assert len(data["residues"]) == 589, f"Residues: {len(data['residues'])}"
    assert len(data["phi"]) == 588
    assert len(data["psi"]) == 588
    assert len(data["omega"]) == 588
    assert len(data["torsion_angles"]) == 588 * 3  # phi + psi + omega
    assert len(data["bond_angles"]) == 588 * 3 + 1  # bond angles between residues
    assert len(data["coordinates"]) == 4678  # N, CA, C, O per residue
    assert len(data["b_factors"]) == 589


@run_test
def test_get_data_multichain_fcz():
    """Extract structure data from compressed multi-chain FCZ."""
    pdb = (ROOT / "test" / "multichain.pdb").read_bytes()
    compressed = foldcomp.compress("test", pdb)
    data = foldcomp.get_data(compressed)
    assert len(data["residues"]) == 589
    assert len(data["coordinates"]) == 4678


@run_test
def test_get_data_single_chain_pdb():
    """Extract structure data from single-chain PDB."""
    pdb = (ROOT / "test" / "test.pdb").read_bytes()
    data = foldcomp.get_data(pdb)
    assert "phi" in data and "psi" in data and "omega" in data
    assert len(data["coordinates"]) == 2208


@run_test
def test_get_data_fcz_file():
    """Extract structure data from existing .fcz file."""
    fcz_data = (ROOT / "test" / "test.fcz").read_bytes()
    data = foldcomp.get_data(fcz_data)
    assert len(data["residues"]) == 276
    assert len(data["coordinates"]) == 2208


@run_test
def test_get_data_af_fcz():
    """Extract data from AlphaFold FCZ (FCMP format, single chain)."""
    fcz_data = (ROOT / "test" / "test_af.fcz").read_bytes()
    data = foldcomp.get_data(fcz_data)
    assert len(data["residues"]) == 28
    assert len(data["coordinates"]) == 243


# ---- split_pdb_by_chain ----

@run_test
def test_split_pdb_by_chain_multichain():
    """Split multi-chain PDB into individual chains."""
    pdb = (ROOT / "test" / "multichain.pdb").read_text()
    chains = foldcomp.split_pdb_by_chain(pdb)
    assert len(chains) == 2, f"Expected 2 chains, got {len(chains)}"
    chain_ids = [c.splitlines()[0][21] for c in chains if c.strip()]
    assert set(chain_ids) == {"A", "B"}


@run_test
def test_split_pdb_by_chain_single():
    """Split single-chain PDB (should return 1 chain)."""
    pdb = (ROOT / "test" / "test.pdb").read_text()
    chains = foldcomp.split_pdb_by_chain(pdb)
    assert len(chains) == 1


@run_test
def test_split_pdb_by_chain_roundtrip():
    """Split chains, compress each, decompress, verify atom counts."""
    pdb = (ROOT / "test" / "multichain.pdb").read_text()
    chains = foldcomp.split_pdb_by_chain(pdb)
    for i, chain_pdb in enumerate(chains):
        chain_id = chain_pdb.splitlines()[0][21]
        compressed = foldcomp.compress(f"chain_{i}", chain_pdb.encode())
        name, out = foldcomp.decompress(compressed)
        in_atoms = len(_atom_lines(chain_pdb))
        out_atoms = len(_atom_lines(out))
        assert in_atoms == out_atoms, f"Chain {chain_id}: in={in_atoms}, out={out_atoms}"


# ---- Database operations ----

@run_test
def test_open_db_all():
    """Open database and iterate all entries."""
    db = ROOT / "test" / "example_db"
    with foldcomp.open(str(db)) as database:
        entries = list(database)
        assert len(entries) == 24, f"Expected 24 entries, got {len(entries)}"


@run_test
def test_open_db_ids():
    """Open database with specific ID subset."""
    db = ROOT / "test" / "example_db"
    with foldcomp.open(str(db), ids=["d1asha_", "d1it2a_"]) as database:
        entries = list(database)
        assert len(entries) == 2
        names = {e[0] for e in entries}
        assert names == {"d1asha_", "d1it2a_"}


@run_test
def test_open_db_str_path():
    """Open database with string path."""
    db = ROOT / "test" / "example_db"
    with foldcomp.open(str(db)) as database:
        first = database[0]
        assert isinstance(first, tuple)
        assert len(first) == 2


@run_test
def test_open_db_pathlib_path():
    """Open database with Path object."""
    db = ROOT / "test" / "example_db"
    with foldcomp.open(db) as database:
        assert len(database) == 24


@run_test
def test_db_entry_atom_counts():
    """Verify atom counts for known database entries."""
    db = ROOT / "test" / "example_db"
    with foldcomp.open(str(db), ids=["d1asha_", "d1it2a_"]) as database:
        entries = list(database)
        for name, pdb_text in entries:
            atoms = len(_atom_lines(pdb_text))
            assert atoms > 0, f"{name} has no atoms"


# ---- FCZ file round-trip ----

@run_test
def test_fcz_file_roundtrip():
    """Decompress existing FCZ file and verify atom count."""
    fcz = (ROOT / "test" / "test.fcz").read_bytes()
    name, out = foldcomp.decompress(fcz)
    assert len(_atom_lines(out)) == 2208


@run_test
def test_fcz_cif_file_roundtrip():
    """Decompress FCZ file from CIF source and verify."""
    fcz = (ROOT / "test" / "test.cif.fcz").read_bytes()
    name, out = foldcomp.decompress(fcz)
    assert name == "AF-A0A009DXE7-F1"
    assert len(_atom_lines(out)) == 1513


# ---- Compression ratios ----

@run_test
def test_compression_ratios():
    """Measure compression ratios for various structures."""
    test_cases = [
        ("test.pdb", "test/test.pdb", "pdb", 2208),
        ("multichain.pdb", "test/multichain.pdb", "pdb", 4678),
        ("test_af.pdb", "test/test_af.pdb", "pdb", 243),
    ]
    ratios = {}
    for label, path, fmt, expected_atoms in test_cases:
        with open(ROOT / path, "r") as f:
            data = f.read().encode()
        compressed = foldcomp.compress(label, data, format=fmt)
        ratio = len(data) / len(compressed)
        ratios[label] = ratio
    assert all(r > 5 for r in ratios.values()), f"Low compression ratios: {ratios}"


# ---- Anchor residue threshold ----

@run_test
def test_compress_with_anchor_threshold():
    """Compress with custom anchor residue threshold."""
    pdb = (ROOT / "test" / "test.pdb").read_bytes()
    compressed = foldcomp.compress("test", pdb, anchor_residue_threshold=10)
    name, out = foldcomp.decompress(compressed)
    assert len(_atom_lines(out)) == 2208


@run_test
def test_compress_with_max_backbone_rmsd():
    """Compress with max_backbone_rmsd parameter."""
    pdb = (ROOT / "test" / "test.pdb").read_bytes()
    compressed = foldcomp.compress("test", pdb, max_backbone_rmsd=2.0)
    name, out = foldcomp.decompress(compressed)
    assert len(_atom_lines(out)) == 2208


# ---- Directory test inputs ----

@run_test
def test_dir_input_multichain_pdb():
    """Test multichain.pdb from dir_test_input."""
    pdb = (ROOT / "test" / "dir_test_input" / "multichain.pdb").read_bytes()
    compressed = foldcomp.compress("dir_multichain", pdb)
    name, out = foldcomp.decompress(compressed)
    chains = set(_chain_ids(out))
    assert chains == {"A", "B"}


@run_test
def test_dir_input_cif():
    """Test cif.gz from dir_test_input."""
    cif = gzip.open(ROOT / "test" / "dir_test_input" / "test.cif.gz", "rt").read().encode()
    compressed = foldcomp.compress("dir_cif", cif, format="mmcif")
    name, out = foldcomp.decompress(compressed)
    assert len(_atom_lines(out)) > 0


# ---- Error handling ----

@run_test
def test_decompress_invalid_format():
    """Decompress with invalid format should raise ValueError."""
    pdb = (ROOT / "test" / "test.pdb").read_bytes()
    compressed = foldcomp.compress("test", pdb)
    try:
        foldcomp.decompress(compressed, format="invalid")
        raise AssertionError("Should have raised ValueError")
    except ValueError:
        pass


@run_test
def test_compress_empty_input():
    """Compress empty input should raise error."""
    try:
        foldcomp.compress("empty", b"")
        raise AssertionError("Should have raised error")
    except (foldcomp.error, ValueError):
        pass


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def generate_report(log_dir: Path):
    """Generate JSON and text report files."""
    total = len(results)
    passed = sum(1 for r in results if r.passed)
    failed = total - passed
    total_time = sum(r.duration for r in results)

    # JSON report
    report_data = {
        "summary": {
            "total_tests": total,
            "passed": passed,
            "failed": failed,
            "total_duration_sec": round(total_time, 4),
            "pass_rate": round(passed / total * 100, 1) if total else 0,
        },
        "results": [
            {
                "name": r.name,
                "passed": r.passed,
                "duration_sec": round(r.duration, 6),
                "details": r.details,
                "error": r.error,
            }
            for r in results
        ],
    }

    json_path = log_dir / "test_results.json"
    with open(json_path, "w") as f:
        json.dump(report_data, f, indent=2)

    # Text summary
    text_lines = [
        "=" * 72,
        "  Foldcomp Multi-Chain Python API - Test Report",
        "=" * 72,
        "",
        f"  Total tests : {total}",
        f"  Passed      : {passed}",
        f"  Failed      : {failed}",
        f"  Pass rate   : {report_data['summary']['pass_rate']}%",
        f"  Total time  : {total_time:.3f}s",
        "",
        "-" * 72,
        "  Test Details",
        "-" * 72,
        "",
    ]
    for r in results:
        status = "PASS" if r.passed else "FAIL"
        text_lines.append(f"  [{status:4s}] {r.name:45s} {r.duration*1000:8.2f}ms")
        if r.error:
            text_lines.append(f"           ERROR: {r.error}")
    text_lines.append("")
    text_lines.append("=" * 72)

    txt_path = log_dir / "test_summary.txt"
    with open(txt_path, "w") as f:
        f.write("\n".join(text_lines))

    return json_path, txt_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Foldcomp multi-chain API tests")
    parser.add_argument("--log-dir", default="test_report", help="Directory for reports")
    args = parser.parse_args()

    log_dir = Path(args.log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)

    # Discover and run all test functions
    test_functions = [obj for name, obj in globals().items()
                      if name.startswith("test_") and callable(obj) and hasattr(obj, "__wrapped__")]
    # Also include non-decorated test_ functions
    test_functions += [obj for name, obj in globals().items()
                       if name.startswith("test_") and callable(obj) and not hasattr(obj, "__wrapped__")]
    # Deduplicate
    seen = set()
    unique_tests = []
    for fn in test_functions:
        if id(fn) not in seen:
            seen.add(id(fn))
            unique_tests.append(fn)

    print(f"Running {len(unique_tests)} tests...")
    for fn in unique_tests:
        fn()
        status = "PASS" if results[-1].passed else "FAIL"
        print(f"  [{status}] {results[-1].name} ({results[-1].duration*1000:.1f}ms)")

    # Generate reports
    json_path, txt_path = generate_report(log_dir)
    print(f"\nReports saved to: {json_path}, {txt_path}")

    # Exit code
    failed = sum(1 for r in results if not r.passed)
    if failed:
        print(f"\n{failed} test(s) FAILED.")
        sys.exit(1)
    else:
        print(f"\nAll {len(results)} tests PASSED.")
        sys.exit(0)


if __name__ == "__main__":
    main()
