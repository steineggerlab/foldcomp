# Foldcomp Multi-Chain Python API — Test Report

> **Generated:** 2026-05-29  
> **Project:** foldcomp (multi-chain support branch)  
> **Python:** 3.13.13 (pytest 9.0.3)  
> **Platform:** macOS arm64 (Darwin 26.0)

---

## 1. Summary

| Metric              | Value     |
|---------------------|-----------|
| Total tests         | **32**    |
| Passed              | **32**    |
| Failed              | **0**     |
| Pass rate           | **100%**  |
| Total execution time| **67.5 ms** |

---

## 2. API Surface Tested

| Function               | Signature                                              | Multi-chain support |
|------------------------|--------------------------------------------------------|---------------------|
| `compress()`           | `compress(name, pdb_content, format, anchor_residue_threshold, max_backbone_rmsd)` | ✅ FCZC container |
| `decompress()`         | `decompress(input, format)`                            | ✅ FCZC container   |
| `get_data()`           | `get_data(input, format)`                              | ✅ All chains       |
| `open()`               | `open(path, ids, decompress, err_on_missing)`          | ✅ Database         |
| `split_pdb_by_chain()` | `split_pdb_by_chain(pdb_str)`                          | ✅ Utility          |

### Format Magics

| Magic | Meaning                          | Multi-chain? |
|-------|----------------------------------|--------------|
| `FCMP` | Legacy single-chain FCZ         | ❌            |
| `FCZC` | New container format            | ✅            |

---

## 3. Test Data Description

### 3.1 PDB Files

| File | Path | Size | Chains | ATOMs | Description |
|------|------|------|--------|-------|-------------|
| `multichain.pdb` | `test/multichain.pdb` | 408KB | A, B | 4,678 | 6PP9 BRAF:MEK1 kinase complex (human) |
| `test.pdb` | `test/test.pdb` | 77KB | A | 2,208 | Single-chain BRAF kinase fragment |
| `test_af.pdb` | `test/test_af.pdb` | 15KB | A | 243 | AlphaFold prediction (A0A0B7P221), pLDDT in B-factors |

#### `test/multichain.pdb` — First 6 ATOM lines
```
ATOM      1  N   SER A 447      24.661 -50.704   4.103  1.00 83.63           N
ATOM      2  CA  SER A 447      25.120 -49.320   4.096  1.00 72.79           C
ATOM      3  C   SER A 447      24.176 -48.426   3.293  1.00 68.83           C
ATOM      4  O   SER A 447      23.167 -48.892   2.762  1.00 69.54           O
ATOM      5  CB  SER A 447      26.544 -49.233   3.541  1.00 73.01           C
ATOM      6  OG  SER A 447      26.704 -50.070   2.409  1.00 65.82           O
```

#### `test/multichain.pdb` — Chain breakdown
| Chain | Chain ID | ATOMs | Residues | Description |
|-------|----------|-------|----------|-------------|
| 1 | A | 2,213 | ~279 | BRAF kinase (human, UniProt P15056) |
| 2 | B | 2,465 | ~310 | MEK1 kinase (human, UniProt Q02750) |

### 3.2 mmCIF Files (gzip compressed)

| File | Path | Chains | ATOMs | Description |
|------|------|--------|-------|-------------|
| `1jwt.cif.gz` | `test/1jwt.cif.gz` | A | 2,429 | Single-chain |
| `1dlp.cif.gz` | `test/1dlp.cif.gz` | A,B,C,D,E,F | 10,364 | 6-chain complex |
| `2ap2.cif.gz` | `test/2ap2.cif.gz` | A,C,P,Q | 3,890 | P-glycoprotein complex |
| `4kuk.cif.gz` | `test/4kuk.cif.gz` | A | 985 | LOV domain (single) |
| `7c2s.cif.gz` | `test/7c2s.cif.gz` | A,B,G,H,I,M | 1,260 | Dengue antibody Fab |
| `1hrn.cif.gz` | `test/1hrn.cif.gz` | A,B | 5,127 | Ribonuclease inhibitor |
| `2gn5.cif.gz` | `test/2gn5.cif.gz` | A | 682 | Small single-chain |
| `1adn.cif.gz` | `test/1adn.cif.gz` | A | 10,234 | DNA complex |
| `1cfg.cif.gz` | `test/1cfg.cif.gz` | A | 3,880 | Single-chain |
| `1lyz.cif.gz` | `test/1lyz.cif.gz` | A | 1,001 | Lysozyme |

### 3.3 FCZ Compressed Files

| File | Path | Size | Magic | Name | ATOMs | Format |
|------|------|------|-------|------|-------|--------|
| `test.fcz` | `test/test.fcz` | 5,223B | FCZC | test (BRAF fragment) | 2,208 | Container |
| `test_af.fcz` | `test/test_af.fcz` | 698B | FCMP | AF prediction A0A0B7P221 | 243 | Legacy |
| `test.cif.fcz` | `test/test.cif.fcz` | 3,113B | FCMP | AF-A0A009DXE7-F1 | 1,513 | Legacy |

### 3.4 Database

| File | Path | Entries | Description |
|------|------|---------|-------------|
| `example_db` | `test/example_db` | 24 | Sample Foldcomp database with `.index`, `.lookup`, `.dbtype`, `.source` files |

#### Sample entries from `example_db`
| ID | ATOMs | Description |
|----|-------|-------------|
| `d1asha_` | 1,240 | SH3 domain |
| `d1it2a_` | 1,189 | Titin fragment |
| `d1b0ba_` | 1,035 | Bacteriorhodopsin |
| `d1cg5a_` | 1,116 | CG5 protein |

### 3.5 Auxiliary Files

| File | Path | Format | Description |
|------|------|--------|-------------|
| `test_af.plddt` | `test/test_af.plddt` | FASTA-like | Single-digit pLDDT scores |
| `test_af.plddt.tsv` | `test/test_af.plddt.tsv` | TSV | Two-digit pLDDT scores |
| `test_fcz.cif` | `test/test_fcz.cif` | mmCIF | Decompressed CIF from test.fcz |
| `test.cif_fcz.cif` | `test/test.cif_fcz.cif` | mmCIF | Decompressed CIF from test.cif.fcz |

### 3.6 Directory Test Inputs

| File | Path | Chains | ATOMs | Format |
|------|------|--------|-------|--------|
| `multichain.pdb` | `test/dir_test_input/multichain.pdb` | A,B | 4,678 | PDB |
| `test.pdb` | `test/dir_test_input/test.pdb` | A | 2,208 | PDB |
| `test.cif.gz` | `test/dir_test_input/test.cif.gz` | A | 1,513 | mmCIF |
| `test_af.pdb` | `test/dir_test_input/test_af.pdb` | A | 243 | PDB |

---

## 4. Test Results by Category

### 4.1 Compress / Decompress (10 tests)

| # | Test | Duration | Status |
|---|------|----------|--------|
| 1 | `test_compress_multichain_pdb` | 1.6ms | ✅ |
| 2 | `test_decompress_multichain_pdb` | 3.4ms | ✅ |
| 3 | `test_decompress_multichain_mmcif` | 5.6ms | ✅ |
| 4 | `test_compress_single_chain_pdb` | 1.4ms | ✅ |
| 5 | `test_compress_af_pdb` | 0.3ms | ✅ |
| 6 | `test_compress_cif` | 2.0ms | ✅ |
| 7 | `test_compress_multichain_cif_1dlp` | 12.5ms | ✅ |
| 8 | `test_compress_multichain_cif_2ap2` | 4.4ms | ✅ |
| 9 | `test_compress_multichain_cif_7c2s` | 2.0ms | ✅ |
| 10 | `test_compress_multichain_cif_1hrn` | 5.4ms | ✅ |

**Key findings:**
- Multi-chain PDB (6PP9) compresses to FCZC container format (magic: `FCZC`)
- 6-chain CIF (1dlp, 10,364 atoms) compresses to 134,876 bytes — all 6 chains preserved
- 4-chain CIF (2ap2) preserves chains A, C, P, Q correctly
- 6-chain CIF (7c2s) preserves chains A, B, G, H, I, M correctly
- mmCIF decompression produces valid CIF with `_atom_site.label_asym_id` chains
- Single-chain PDBs produce FCZC containers when > ~25 residues

### 4.2 get_data — Structure Data Extraction (5 tests)

| # | Test | Duration | Status |
|---|------|----------|--------|
| 11 | `test_get_data_multichain_pdb` | 1.5ms | ✅ |
| 12 | `test_get_data_multichain_fcz` | 2.5ms | ✅ |
| 13 | `test_get_data_single_chain_pdb` | 0.6ms | ✅ |
| 14 | `test_get_data_fcz_file` | 0.6ms | ✅ |
| 15 | `test_get_data_af_fcz` | 0.3ms | ✅ |

**Multi-chain get_data output (6PP9 BRAF:MEK1):**

| Key | Length | Description |
|-----|--------|-------------|
| `phi` | 588 | Phi torsion angles (C-N-CA-C) |
| `psi` | 588 | Psi torsion angles (N-CA-C-N) |
| `omega` | 588 | Omega torsion angles (CA-C-N-CA) |
| `torsion_angles` | 1,764 | Concatenated phi+psi+omega (588×3) |
| `bond_angles` | 1,765 | Backbone bond angles (CA-C-N, N-CA-C) |
| `residues` | 589 | Amino acid sequence string |
| `coordinates` | 4,678 | 3D coordinates (N, CA, C, O per residue) |
| `b_factors` | 589 | B-factors / pLDDT scores |

**Key findings:**
- `get_data()` concatenates data across all chains in the structure
- Residue count (589) = Chain A (~279) + Chain B (~310)
- Coordinates (4,678) = 589 residues × ~8 atoms/residue backbone

### 4.3 split_pdb_by_chain (3 tests)

| # | Test | Duration | Status |
|---|------|----------|--------|
| 16 | `test_split_pdb_by_chain_multichain` | 0.8ms | ✅ |
| 17 | `test_split_pdb_by_chain_single` | 0.3ms | ✅ |
| 18 | `test_split_pdb_by_chain_roundtrip` | 3.4ms | ✅ |

**Key findings:**
- `split_pdb_by_chain("6PP9")` → 2 chains: A (2,213 atoms), B (2,465 atoms)
- Each split chain compresses and decompresses independently with 100% atom fidelity

### 4.4 Database Operations (5 tests)

| # | Test | Duration | Status |
|---|------|----------|--------|
| 19 | `test_open_db_all` | 6.5ms | ✅ |
| 20 | `test_open_db_ids` | 0.7ms | ✅ |
| 21 | `test_open_db_str_path` | 0.3ms | ✅ |
| 22 | `test_open_db_pathlib_path` | 0.03ms | ✅ |
| 23 | `test_db_entry_atom_counts` | 0.7ms | ✅ |

**Key findings:**
- Database with 24 entries opens and iterates correctly
- ID-based filtering (`ids=["d1asha_", "d1it2a_"]`) returns exactly 2 entries
- Both string paths and `pathlib.Path` objects work as `path` argument
- Entries decompress to valid PDB with correct atom counts

### 4.5 FCZ File Round-Trip (2 tests)

| # | Test | Duration | Status |
|---|------|----------|--------|
| 24 | `test_fcz_file_roundtrip` | 0.7ms | ✅ |
| 25 | `test_fcz_cif_file_roundtrip` | 0.7ms | ✅ |

**Key findings:**
- `test.fcz` (FCZC container, 5,223B) → decompresses to 2,208 atoms
- `test.cif.fcz` (FCMP legacy, 3,113B) → name: "AF-A0A009DXE7-F1", 1,513 atoms

### 4.6 Compression Ratios (1 test)

| # | Test | Duration | Status |
|---|------|----------|--------|
| 26 | `test_compression_ratios` | 1.9ms | ✅ |

**Measured compression ratios:**

| File | Input size | Compressed | Ratio |
|------|-----------|------------|-------|
| `test.pdb` | 77KB | 5,034B | **15.3x** |
| `multichain.pdb` | 408KB | 29,070B | **14.0x** |
| `test_af.pdb` | 15KB | 638B | **23.7x** |

### 4.7 Advanced Parameters (2 tests)

| # | Test | Duration | Status |
|---|------|----------|--------|
| 27 | `test_compress_with_anchor_threshold` | 1.2ms | ✅ |
| 28 | `test_compress_with_max_backbone_rmsd` | 1.3ms | ✅ |

### 4.8 Directory Test Inputs (2 tests)

| # | Test | Duration | Status |
|---|------|----------|--------|
| 29 | `test_dir_input_multichain_pdb` | 2.8ms | ✅ |
| 30 | `test_dir_input_cif` | 1.8ms | ✅ |

### 4.9 Error Handling (2 tests)

| # | Test | Duration | Status |
|---|------|----------|--------|
| 31 | `test_decompress_invalid_format` | 0.5ms | ✅ |
| 32 | `test_compress_empty_input` | 0.02ms | ✅ |

---

## 5. Multi-Chain Support Architecture

### 5.1 Container Format (FCZC)

When a multi-chain structure is compressed, foldcomp uses a new **container format**:

```
FCZC [CompressedFileHeader] [nFragments: uint32] [title] [fragments...]
```

Each fragment contains:
- `kind` (uint8): `0`=FCZ encoded, `2`=raw atoms
- `model` (int16): model number
- `chain` (string): chain identifier
- `payload` (bytes): compressed or raw data

### 5.2 Raw Atom Fallback

Regions that cannot be encoded with foldcomp (non-standard residues, insertion codes, OXT atoms beyond C-terminus) fall back to raw atom storage within the container.

### 5.3 Chain Identification

The `identifyChains()` function in `structure_codec.cpp` groups atoms by chain ID, then each chain is independently processed through:
1. `identifyDiscontinousResInd()` — find continuous residue spans
2. `identifyBackboneRegions()` — find backbone-encodable regions
3. `regionNeedsRawFallback()` — check if region needs raw storage
4. Compression or raw serialization

---

## 6. Files in This Test Report

| File | Purpose |
|------|---------|
| `TEST_REPORT.md` | This report |
| `test_multichain_api.py` | Comprehensive test suite (32 tests) |
| `test_results.json` | Structured JSON results |
| `test_summary.txt` | Human-readable text summary |
| `test_data/` | Test data directory (symlinked to `../test/`) |

---

## 7. How to Re-run

```bash
cd /Users/hyunbin/projects/lab/foldcomp/maintenance/multichain_support/foldcomp

# Using the virtual environment Python
.venv/bin/python3 test_report/test_multichain_api.py --log-dir test_report

# Or with pytest
.venv/bin/pytest test/test_foldcomp.py -v
```

---

## 8. Conclusions

✅ **All 32 tests pass** — The multi-chain support implementation is working correctly across all tested scenarios:

1. **PDB multi-chain** (2-chain 6PP9): Compress → FCZC container → Decompress preserves both chains A and B
2. **CIF multi-chain** (2, 4, and 6 chains): All chains correctly identified and preserved through round-trip
3. **get_data()**: Extracts torsion angles, coordinates, residues across all chains
4. **split_pdb_by_chain()**: Correctly splits multi-chain PDBs for per-chain processing
5. **Database operations**: Full iteration and ID-based subsetting work correctly
6. **Compression ratios**: 14-24x compression achieved, comparable to single-chain performance
7. **Error handling**: Invalid format and empty input raise appropriate exceptions
