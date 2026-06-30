#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 1 ]; then
  echo "usage: $0 <foldcomp-binary>" >&2
  exit 2
fi

FOLDCOMP_BIN=$1
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TEST_DIR=$SCRIPT_DIR
TMP_ROOT=${FOLDCOMP_SMOKE_TMPDIR:-/tmp}
mkdir -p "$TMP_ROOT"

assert_close() {
  local observed=$1
  local expected=$2
  local tolerance=$3
  awk -v check="$observed" -v target="$expected" -v tol="$tolerance" \
    'BEGIN { diff = check - target; if (diff < 0) diff = -diff; if (diff > tol) { print check "!=" target; exit 1 } }'
}

assert_rmsd_exact() {
  local source=$1
  local roundtrip=$2
  local expected_bb=$3
  local expected_aa=$4
  local tolerance=${5:-0.001}
  local line bb aa
  line=$("$FOLDCOMP_BIN" rmsd "$source" "$roundtrip")
  bb=$(printf '%s\n' "$line" | cut -f5)
  aa=$(printf '%s\n' "$line" | cut -f6)
  assert_close "$bb" "$expected_bb" "$tolerance"
  assert_close "$aa" "$expected_aa" "$tolerance"
}

assert_rmsd_at_most() {
  local source=$1
  local roundtrip=$2
  local max_bb=$3
  local max_aa=$4
  local line bb aa
  line=$("$FOLDCOMP_BIN" rmsd "$source" "$roundtrip")
  bb=$(printf '%s\n' "$line" | cut -f5)
  aa=$(printf '%s\n' "$line" | cut -f6)
  awk -v bb="$bb" -v max_bb="$max_bb" -v aa="$aa" -v max_aa="$max_aa" \
    'BEGIN {
      if (bb > max_bb) { print "backbone " bb ">" max_bb; exit 1 }
      if (aa > max_aa) { print "all-atom " aa ">" max_aa; exit 1 }
    }'
}

run_reference_cases() {
  rm -f "$TEST_DIR/test.fcz" "$TEST_DIR/test_fcz.cif" "$TEST_DIR/test.cif.fcz" "$TEST_DIR/test.cif_fcz.cif"

  "$FOLDCOMP_BIN" compress "$TEST_DIR/test.pdb"
  "$FOLDCOMP_BIN" decompress "$TEST_DIR/test.fcz" "$TEST_DIR/test_fcz.cif"
  assert_rmsd_exact "$TEST_DIR/test.pdb" "$TEST_DIR/test_fcz.cif" "0.0453689" "0.0833658"

  "$FOLDCOMP_BIN" compress "$TEST_DIR/test.cif.gz"
  "$FOLDCOMP_BIN" decompress -a "$TEST_DIR/test.cif.fcz" "$TEST_DIR/test.cif_fcz.cif"
  assert_rmsd_exact "$TEST_DIR/test.cif.gz" "$TEST_DIR/test.cif_fcz.cif" "0.0569663" "0.128881"
}

run_regression_case() {
  local name=$1
  local max_default_bb=$2
  local max_default_aa=$3
  local max_spill_bb=$4
  local max_spill_aa=$5
  local source="$TEST_DIR/${name}.cif.gz"
  local default_base="${TMP_ROOT}/foldcomp-smoke-${name}-default"
  local max_base="${TMP_ROOT}/foldcomp-smoke-${name}-max"

  rm -f "${default_base}" "${default_base}." "${default_base}.cif"
  "$FOLDCOMP_BIN" compress "$source" "$default_base"
  "$FOLDCOMP_BIN" decompress "${default_base}." "${default_base}.cif"
  assert_rmsd_at_most "$source" "${default_base}.cif" "$max_default_bb" "$max_default_aa"

  rm -f "${max_base}" "${max_base}." "${max_base}.cif"
  "$FOLDCOMP_BIN" compress --max-backbone-rmsd 0.2 "$source" "$max_base"
  "$FOLDCOMP_BIN" decompress "${max_base}." "${max_base}.cif"
  assert_rmsd_at_most "$source" "${max_base}.cif" "$max_spill_bb" "$max_spill_aa"
}

run_reference_cases
run_regression_case "1cfg" "0.2" "0.7" "0.2" "0.7"
run_regression_case "1adn" "0.35" "0.4" "0.2" "0.5"
run_regression_case "1dlp" "0.2" "0.6" "0.2" "0.5"
run_regression_case "1hrn" "0.12" "0.2" "0.2" "0.5"
run_regression_case "1jwt" "0.16" "0.35" "0.2" "0.5"
run_regression_case "1lyz" "0.55" "0.75" "0.2" "0.5"
run_regression_case "2ap2" "0.15" "0.2" "0.2" "0.5"
run_regression_case "2gn5" "1.5" "1.8" "0.2" "0.8"
run_regression_case "4kuk" "0.1" "0.15" "0.2" "0.5"
run_regression_case "7c2s" "0.01" "0.01" "0.2" "0.5"
