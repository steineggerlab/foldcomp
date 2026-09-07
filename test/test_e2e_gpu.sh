#!/bin/bash
# test/test_e2e_gpu.sh
# End-to-end GPU decompression test suite.
#
# Tests the compress→decompress roundtrip accuracy (RMSD), batch decompression
# from a database, and correctness across backbone reconstruction algorithms.
#
# Usage:
#   ./test/test_e2e_gpu.sh [binary_path]
#   binary_path defaults to ./build/foldcomp
#
# Exit code: 0 on all tests passing, non-zero on failure.

set -e

FC="${1:-./build_gpu/foldcomp}"
TMPDIR_BASE="${TMPDIR:-/tmp}/foldcomp_e2e_$$"
PASS=0
FAIL=0

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

cleanup() {
  rm -rf "$TMPDIR_BASE"
}
trap cleanup EXIT

die() {
  echo "[FAIL] $*" >&2
  FAIL=$((FAIL + 1))
}

pass() {
  echo "[PASS] $*"
  PASS=$((PASS + 1))
}

# Check that RMSD (field 5 of foldcomp rmsd output) is at or below a ceiling.
# Unlike a fixed-golden RMSD check, this takes no expected value: it asserts a relationship
# (e.g. "GPU output tracks the CPU reference to within X Angstrom") rather
# than pinning an exact, arch/fp32-rounding-sensitive constant.
# Usage: check_rmsd_max <ref_pdb> <query_pdb> <max> <label>
check_rmsd_max() {
  local ref="$1" query="$2" max="$3" label="$4"
  local rmsd
  rmsd=$("$FC" rmsd "$ref" "$query" 2>&1 | cut -f5)
  if [ -z "$rmsd" ]; then
    die "$label: rmsd command failed (empty output)"
    return
  fi
  awk -v got="$rmsd" -v max="$max" -v lbl="$label" '
        BEGIN {
            if (got > max) {
                printf "[FAIL] %s: RMSD %.6f exceeds max %.4f\n", lbl, got, max > "/dev/stderr"
                exit 1
            }
        }
    ' && pass "$label: RMSD $rmsd (max $max)" ||
    {
      FAIL=$((FAIL + 1))
      PASS=$((PASS > 0 ? PASS - 1 : 0))
    }
}

# Check that file exists and has non-zero size.
check_file() {
  local f="$1" label="$2"
  if [ -s "$f" ]; then
    pass "$label: file exists ($f)"
  else
    die "$label: file missing or empty ($f)"
  fi
}

# Count ATOM/HETATM lines in a PDB file.
atom_count() {
  grep -c "^ATOM\|^HETATM" "$1" 2>/dev/null || echo 0
}

check_atom_count() {
  local pdb="$1" expected="$2" label="$3"
  local got
  got=$(atom_count "$pdb")
  if [ "$got" -eq "$expected" ]; then
    pass "$label: atom count $got"
  else
    die "$label: atom count $got, expected $expected"
  fi
}

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TEST_DIR="$REPO_ROOT/test"

mkdir -p "$TMPDIR_BASE"

if [ ! -x "$FC" ]; then
  echo "ERROR: binary not found: $FC" >&2
  exit 1
fi

# This suite's whole value is asserting "GPU matches CPU" — against a CPU-only
# (-DBUILD_CUDA=OFF) binary, --gpu silently falls back to --cpu, so every such
# comparison degrades to CPU-vs-CPU and the suite goes green while testing
# nothing (see MR !1 note_67749589). Guard on the build flavour foldcomp
# --version reports instead of running anything against the wrong binary.
if ! "$FC" --version | grep -q '(cuda)'; then
  echo "ERROR: $FC is not a CUDA-enabled build (rebuild with -DBUILD_CUDA=ON);" \
       "GPU-vs-CPU comparisons would silently be CPU-vs-CPU." >&2
  exit 1
fi

echo "=== foldcomp GPU end-to-end tests ==="
echo "Binary: $FC"
echo "Temp:   $TMPDIR_BASE"
echo ""

# ---------------------------------------------------------------------------
# Test 1: PDB roundtrip — compress then batch-decompress (single protein)
# ---------------------------------------------------------------------------
echo "--- Test 1: PDB compress→decompress roundtrip ---"
{
  FCZ="$TMPDIR_BASE/test.fcz"
  IN_DIR="$TMPDIR_BASE/t1_in"
  OUT_DIR="$TMPDIR_BASE/t1_out"
  mkdir -p "$IN_DIR"

  "$FC" compress -y "$TEST_DIR/test.pdb" "$FCZ"
  check_file "$FCZ" "T1 compress"

  cp "$FCZ" "$IN_DIR/"
  "$FC" decompress --write-to-disk -y -t 1 "$IN_DIR" "$OUT_DIR"
  DEC="$OUT_DIR/test.pdb"
  check_file "$DEC" "T1 decompress output"

  # test.pdb: 2208 atoms. The GPU path appends 1 OXT slot (may be null if no OXT
  # in source); RMSD uses only non-null atoms so 2208-atom RMSD is valid.
  check_atom_count "$DEC" 2208 "T1 atom count"

  CPU_DEC="$OUT_DIR/test_cpu.pdb"
  "$FC" decompress --cpu -y "$FCZ" "$CPU_DEC"
  check_file "$CPU_DEC" "T1 CPU decompress output"
  # Sanity bound on the reference (CPU) path, not a GPU-specific golden.
  check_rmsd_max "$TEST_DIR/test.pdb" "$CPU_DEC" 0.2 "T1 CPU roundtrip sanity"
  # The property that matters: GPU reproduces the CPU reference implementation.
  check_rmsd_max "$CPU_DEC" "$DEC" 0.01 "T1 GPU matches CPU"
}

# ---------------------------------------------------------------------------
# Test 2: CIF roundtrip — compress gzipped CIF then batch-decompress
# ---------------------------------------------------------------------------
echo ""
echo "--- Test 2: gzipped CIF compress→decompress roundtrip ---"
{
  FCZ="$TMPDIR_BASE/test_cif.fcz"
  IN_DIR="$TMPDIR_BASE/t2_in"
  OUT_DIR="$TMPDIR_BASE/t2_out"
  mkdir -p "$IN_DIR"

  "$FC" compress -y "$TEST_DIR/test.cif.gz" "$FCZ"
  check_file "$FCZ" "T2 compress"

  cp "$FCZ" "$IN_DIR/"
  "$FC" decompress --write-to-disk -a -y -t 1 "$IN_DIR" "$OUT_DIR"
  DEC="$OUT_DIR/test_cif.pdb"
  check_file "$DEC" "T2 decompress output"
  check_atom_count "$DEC" 1513 "T2 atom count"

  CPU_DEC="$OUT_DIR/test_cif_cpu.pdb"
  "$FC" decompress --cpu -a -y "$FCZ" "$CPU_DEC"
  check_file "$CPU_DEC" "T2 CPU decompress output"
  check_rmsd_max "$TEST_DIR/test.cif.gz" "$CPU_DEC" 0.2 "T2 CPU roundtrip sanity"
  check_rmsd_max "$CPU_DEC" "$DEC" 0.01 "T2 GPU matches CPU"
}

# ---------------------------------------------------------------------------
# Test 3: Database decompress — verify all 24 entries are produced, and that
# every one of them matches the CPU reference (atom count + RMSD). A plain
# file/structure count (as this test used to do) passes even when individual
# entries are corrupted — e.g. a per-residue atom-count bug that overruns
# into the next structure's output in the same batch — since the total file
# count is unaffected. See MR !1 note_67749575.
# ---------------------------------------------------------------------------
echo ""
echo "--- Test 3: Database decompress (example_db) ---"
{
  OUT_CPU="$TMPDIR_BASE/t3_cpu"
  OUT_GPU="$TMPDIR_BASE/t3_gpu"
  "$FC" decompress --cpu --write-to-disk -y "$TEST_DIR/example_db" "$OUT_CPU"
  "$FC" decompress --gpu --write-to-disk -y "$TEST_DIR/example_db" "$OUT_GPU"

  EXPECTED_COUNT=24
  ACTUAL_COUNT=$(ls "$OUT_GPU"/*.pdb 2>/dev/null | wc -l)
  if [ "$ACTUAL_COUNT" -eq "$EXPECTED_COUNT" ]; then
    pass "T3 DB decompress: $ACTUAL_COUNT structures produced"
  else
    die "T3 DB decompress: $ACTUAL_COUNT structures, expected $EXPECTED_COUNT"
  fi

  # Spot-check a few specific structures exist and are non-empty
  for name in d1asha_ d1it2a_ d1hlba_; do
    check_file "$OUT_GPU/${name}.pdb" "T3 $name"
  done

  # Per-entry: GPU must match CPU on both atom count and coordinates. This is
  # the only test covering a multi-structure batch, which is precisely where
  # cross-structure corruption from a per-residue atom-count bug shows up.
  for cpu_pdb in "$OUT_CPU"/*.pdb; do
    b=$(basename "$cpu_pdb")
    gpu_pdb="$OUT_GPU/$b"
    cpu_atoms=$(atom_count "$cpu_pdb")
    gpu_atoms=$(atom_count "$gpu_pdb")
    if [ "$cpu_atoms" -eq "$gpu_atoms" ]; then
      pass "T3 $b atom count matches CPU: $cpu_atoms"
    else
      die "T3 $b atom count differs: CPU $cpu_atoms, GPU $gpu_atoms"
      continue # RMSD comparison needs matching atom counts to be meaningful
    fi
    check_rmsd_max "$cpu_pdb" "$gpu_pdb" 0.01 "T3 $b GPU matches CPU"
  done
}

# ---------------------------------------------------------------------------
# Test 4: Backbone reconstruction correctness (segmented algorithm, default)
# ---------------------------------------------------------------------------
echo ""
echo "--- Test 4: Backbone algorithm correctness ---"
{
  FCZ="$TMPDIR_BASE/test_algo.fcz"
  "$FC" compress -y "$TEST_DIR/test.pdb" "$FCZ"

  IN_DIR="$TMPDIR_BASE/t4_in_segmented"
  OUT_DIR="$TMPDIR_BASE/t4_out_segmented"
  mkdir -p "$IN_DIR"
  cp "$FCZ" "$IN_DIR/"
  "$FC" decompress --write-to-disk -y -t 1 "$IN_DIR" "$OUT_DIR"
  DEC="$OUT_DIR/test_algo.pdb"
  check_file "$DEC" "T4 segmented output"

  CPU_DEC="$OUT_DIR/test_algo_cpu.pdb"
  "$FC" decompress --cpu -y "$FCZ" "$CPU_DEC"
  check_file "$CPU_DEC" "T4 CPU decompress output"
  check_rmsd_max "$TEST_DIR/test.pdb" "$CPU_DEC" 0.2 "T4 CPU roundtrip sanity"
  check_rmsd_max "$CPU_DEC" "$DEC" 0.01 "T4 segmented GPU matches CPU"
}

# ---------------------------------------------------------------------------
# Test 5: Pre-existing FCZ decompression (repo test files)
# ---------------------------------------------------------------------------
echo ""
echo "--- Test 5: Decompress pre-existing FCZ (test_af.fcz) ---"
{
  IN_DIR="$TMPDIR_BASE/t5_in"
  OUT_DIR="$TMPDIR_BASE/t5_out"
  mkdir -p "$IN_DIR"
  cp "$TEST_DIR/test_af.fcz" "$IN_DIR/"
  "$FC" decompress --write-to-disk -y -t 1 "$IN_DIR" "$OUT_DIR"
  check_file "$OUT_DIR/test_af.pdb" "T5 decompress output"
  # 243 atoms in original. GPU path may append OXT/null atoms; just check >= 243
  COUNT=$(atom_count "$OUT_DIR/test_af.pdb")
  if [ "$COUNT" -ge 243 ]; then
    pass "T5 atom count: $COUNT (>= 243)"
  else
    die "T5 atom count: $COUNT (< 243)"
  fi
}

# ---------------------------------------------------------------------------
# Test 6: check --no-fused-disc produces same RMSD (unfused discretizer path)
# ---------------------------------------------------------------------------
echo ""
echo "--- Test 6: Unfused discretizer kernel ---"
{
  FCZ="$TMPDIR_BASE/test_disc.fcz"
  "$FC" compress -y "$TEST_DIR/test.pdb" "$FCZ"
  IN_DIR="$TMPDIR_BASE/t6_in"
  OUT_DIR="$TMPDIR_BASE/t6_out"
  mkdir -p "$IN_DIR"
  cp "$FCZ" "$IN_DIR/"
  "$FC" decompress --write-to-disk --no-fused-disc -y -t 1 "$IN_DIR" "$OUT_DIR"
  DEC="$OUT_DIR/test_disc.pdb"
  check_file "$DEC" "T6 output"

  CPU_DEC="$OUT_DIR/test_disc_cpu.pdb"
  "$FC" decompress --cpu -y "$FCZ" "$CPU_DEC"
  check_file "$CPU_DEC" "T6 CPU decompress output"
  check_rmsd_max "$TEST_DIR/test.pdb" "$CPU_DEC" 0.2 "T6 CPU roundtrip sanity"
  check_rmsd_max "$CPU_DEC" "$DEC" 0.01 "T6 no-fused-disc GPU matches CPU"
}

# ---------------------------------------------------------------------------
# Test 7: --no-mmap path produces identical output to mmap path
# ---------------------------------------------------------------------------
echo ""
echo "--- Test 7: --no-mmap directory decompression ---"
{
  FCZ="$TMPDIR_BASE/test_nommap.fcz"
  "$FC" compress -y "$TEST_DIR/test.pdb" "$FCZ"

  IN_DIR="$TMPDIR_BASE/t7_in"
  OUT_MMAP="$TMPDIR_BASE/t7_out_mmap"
  OUT_FREAD="$TMPDIR_BASE/t7_out_fread"
  mkdir -p "$IN_DIR"
  cp "$FCZ" "$IN_DIR/"

  "$FC" decompress --write-to-disk -y -t 1 "$IN_DIR" "$OUT_MMAP"
  "$FC" decompress --write-to-disk --no-mmap -y -t 1 "$IN_DIR" "$OUT_FREAD"

  check_file "$OUT_FREAD/test_nommap.pdb" "T7 no-mmap output exists"
  check_atom_count "$OUT_FREAD/test_nommap.pdb" 2208 "T7 no-mmap atom count"

  CPU_DEC="$TMPDIR_BASE/t7_cpu.pdb"
  "$FC" decompress --cpu -y "$FCZ" "$CPU_DEC"
  check_file "$CPU_DEC" "T7 CPU decompress output"
  check_rmsd_max "$TEST_DIR/test.pdb" "$CPU_DEC" 0.2 "T7 CPU roundtrip sanity"
  check_rmsd_max "$CPU_DEC" "$OUT_FREAD/test_nommap.pdb" 0.01 "T7 no-mmap GPU matches CPU"

  # Byte-identical with mmap output
  if cmp -s "$OUT_MMAP/test_nommap.pdb" "$OUT_FREAD/test_nommap.pdb"; then
    pass "T7 no-mmap output byte-identical to mmap output"
  else
    die "T7 no-mmap output differs from mmap output"
  fi
}

# ---------------------------------------------------------------------------
# Test 8: FCZC multi-chain fragments cross GPU batches without losing metadata
# ---------------------------------------------------------------------------
echo ""
echo "--- Test 8: FCZC multi-chain GPU assembly across batches ---"
{
  FCZ="$TMPDIR_BASE/multichain.fcz"
  IN_DIR="$TMPDIR_BASE/t8_in"
  OUT_DIR="$TMPDIR_BASE/t8_out"
  mkdir -p "$IN_DIR"

  "$FC" compress -y "$TEST_DIR/multichain.pdb" "$FCZ"
  if [ "$(head -c 4 "$FCZ")" = "FCZC" ]; then
    pass "T8 compression produced FCZC container"
  else
    die "T8 compression did not produce FCZC container"
  fi

  cp "$FCZ" "$IN_DIR/"
  # One fragment per GPU batch ensures logical-container assembly spans batches.
  TIMING_OUTPUT=$("$FC" decompress --write-to-disk --max-batch-size 1 --time -y -t 1 "$IN_DIR" "$OUT_DIR")
  if printf '%s\n' "$TIMING_OUTPUT" | grep -q $'multichain.fcz\t'; then
    die "T8 container used the whole-input CPU fallback"
  else
    pass "T8 container bypassed the whole-input CPU fallback"
  fi
  DEC="$OUT_DIR/multichain.pdb"
  check_file "$DEC" "T8 decompression output"

  EXPECTED_ATOMS=$(atom_count "$TEST_DIR/multichain.pdb")
  check_atom_count "$DEC" "$EXPECTED_ATOMS" "T8 atom count"
  CHAINS=$(awk '/^(ATOM|HETATM)/ {print substr($0,22,1)}' "$DEC" | sort -u | tr '\n' ' ')
  if [ "$CHAINS" = "A B " ]; then
    pass "T8 chain identifiers preserved"
  else
    die "T8 chain identifiers: '$CHAINS', expected 'A B '"
  fi

  CPU_DEC="$TMPDIR_BASE/t8_cpu.pdb"
  "$FC" decompress --cpu -y "$FCZ" "$CPU_DEC"
  IDENTITY_CHECK=$("$FC" rmsd "$CPU_DEC" "$DEC" 2>&1)
  if printf '%s\n' "$IDENTITY_CHECK" | grep -q "identity mismatch"; then
    die "T8 GPU residue/atom identities differ from CPU output"
  else
    pass "T8 GPU residue/atom identities match CPU output"
  fi
}

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
echo "============================="
echo "Results: $PASS passed, $FAIL failed"
echo "============================="

if [ "$FAIL" -gt 0 ]; then
  exit 1
fi
exit 0
