# Test Data Description

## Directory Layout

```
test/
├── multichain.pdb          # Multi-chain PDB (6PP9 BRAF:MEK1)
├── test.pdb                # Single-chain PDB (BRAF fragment)
├── test_af.pdb             # AlphaFold predicted PDB
├── test_af.plddt           # FASTA-like pLDDT scores
├── test_af.plddt.tsv       # TSV pLDDT scores
├── test.fcz                # Compressed FCZ (container format)
├── test_af.fcz             # Compressed FCZ (legacy format)
├── test.cif.fcz            # Compressed FCZ (from CIF source)
├── test_fcz.cif            # mmCIF from test.fcz
├── test.cif_fcz.cif        # mmCIF from test.cif.fcz
├── test.cif.gz             # mmCIF gzip
├── 1jwt.cif.gz             # mmCIF (single-chain)
├── 1dlp.cif.gz             # mmCIF (6 chains: A,B,C,D,E,F)
├── 2ap2.cif.gz             # mmCIF (4 chains: A,C,P,Q)
├── 4kuk.cif.gz             # mmCIF (single-chain)
├── 7c2s.cif.gz             # mmCIF (6 chains: A,B,G,H,I,M)
├── 1hrn.cif.gz             # mmCIF (2 chains: A,B)
├── 2gn5.cif.gz             # mmCIF (single-chain)
├── 1adn.cif.gz             # mmCIF (single-chain)
├── 1cfg.cif.gz             # mmCIF (single-chain)
├── 1lyz.cif.gz             # mmCIF (single-chain)
├── example_db              # Foldcomp database (24 entries)
│   ├── example_db          # Main database file
│   ├── example_db.index    # Index file
│   ├── example_db.lookup   # ID lookup file
│   ├── example_db.dbtype   # Database type
│   └── example_db.source   # Source metadata
├── example_db.subset       # Subset database file
├── pdb_list.txt            # PDB file list
├── fcz_list.txt            # FCZ file list
├── dir_test_input/         # Directory of test inputs
│   ├── multichain.pdb
│   ├── test.pdb
│   ├── test.cif.gz
│   └── test_af.pdb
├── tar_test_input.tar      # Tar archive of test inputs
├── gz_test_input.tar.gz    # Gzipped tar archive
└── test_foldcomp.py        # Existing pytest tests
```

---

## Detailed File Descriptions

### PDB Files

#### `test/multichain.pdb` (408KB)

**Description:** Crystal structure of BRAF:MEK1 complex (PDB ID: 6PP9)
- **Resolution:** 2.59 Å
- **Method:** X-ray diffraction
- **Chains:** 2 (A = BRAF, B = MEK1)

| Chain | Protein | UniProt | Residues | ATOMs |
|-------|---------|---------|----------|-------|
| A | BRAF kinase (human) | P15056 (BRAF_HUMAN) | ~279 | 2,213 |
| B | MEK1 kinase (human) | Q02750 (MP2K1_HUMAN) | ~310 | 2,465 |

**Ligands:** ANP (phosphoaminophosphonic acid-adenylate ester), Mg²⁺, Cl⁻, LCJ (drug inhibitor), SO₄²⁻, GOL (glycerol)

**First lines:**
```
HEADER    TRANSFERASE                             05-JUL-19   6PP9
TITLE     CRYSTAL STRUCTURE OF BRAF:MEK1 COMPLEX
...
ATOM      1  N   SER A 447      24.661 -50.704   4.103  1.00 83.63           N
ATOM      2  CA  SER A 447      25.120 -49.320   4.096  1.00 72.79           C
ATOM      3  C   SER A 447      24.176 -48.426   3.293  1.00 68.83           C
```

#### `test/test.pdb` (77KB)

**Description:** Single-chain BRAF kinase fragment (Chain A only from 6PP9)
- **Chain:** A
- **ATOMs:** 2,208
- **Residues:** 276

#### `test/test_af.pdb` (15KB)

**Description:** AlphaFold v2.0 prediction for uncharacterized protein RDT1 (A0A0B7P221)
- **Chain:** A
- **ATOMs:** 243
- **Residues:** 28
- **B-factors:** pLDDT scores (range: 56.9 – 81.4)

---

### FCZ Compressed Files

#### `test/test.fcz` (5,223 bytes)

- **Magic:** `FCZC` (container format)
- **Name:** test (BRAF fragment)
- **ATOMs:** 2,208
- **Chains:** A
- **Format:** Container (multi-chain capable, single-chain content)

#### `test/test_af.fcz` (698 bytes)

- **Magic:** `FCMP` (legacy format)
- **Name:** ALPHAFOLD V2.0 PREDICTION...
- **ATOMs:** 243
- **Chains:** A
- **Format:** Legacy single-chain FCZ

#### `test/test.cif.fcz` (3,113 bytes)

- **Magic:** `FCMP` (legacy format)
- **Name:** AF-A0A009DXE7-F1
- **ATOMs:** 1,513
- **Chains:** A
- **Format:** Legacy single-chain FCZ

---

### mmCIF Files (gzip compressed)

| File | Chains | ATOMs | Description |
|------|--------|-------|-------------|
| `1jwt.cif.gz` | A | 2,429 | Single-chain |
| `1dlp.cif.gz` | A,B,C,D,E,F | 10,364 | 6-chain complex |
| `2ap2.cif.gz` | A,C,P,Q | 3,890 | P-glycoprotein complex |
| `4kuk.cif.gz` | A | 985 | LOV photoreceptor domain |
| `7c2s.cif.gz` | A,B,G,H,I,M | 1,260 | Dengue serotype antibody Fab |
| `1hrn.cif.gz` | A,B | 5,127 | Ribonuclease inhibitor (dimer) |
| `2gn5.cif.gz` | A | 682 | Small single-chain |
| `1adn.cif.gz` | A | 10,234 | DNA complex |
| `1cfg.cif.gz` | A | 3,880 | Single-chain |
| `1lyz.cif.gz` | A | 1,001 | Lysozyme |

---

### Database

#### `test/example_db` (24 entries)

**Files:**
- `example_db` — Main database binary
- `example_db.index` — Entry index
- `example_db.lookup` — ID-to-index lookup
- `example_db.dbtype` — Database type marker
- `example_db.source` — Source metadata
- `example_db.subset` — Subset database (for createsubdb tests)

**Sample entries:**

| ID | ATOMs | Description |
|----|-------|-------------|
| d1asha_ | 1,240 | SH3 domain |
| d1it2a_ | 1,189 | Titin fragment |
| d1b0ba_ | 1,035 | Bacteriorhodopsin |
| d1cg5a_ | 1,116 | CG5 protein |

---

### Auxiliary Files

#### `test/test_af.plddt` (47 bytes)
FASTA-like format with single-digit pLDDT confidence scores:
```
>test/test_af.fcz
5666777787777787777777776665
```

#### `test/test_af.plddt.tsv` (TSV)
Tab-separated file with pLDDT scores per residue (2-digit precision).

#### `test/test_fcz.cif`
mmCIF output from decompressing `test.fcz` with format="mmcif".

#### `test/test.cif_fcz.cif`
mmCIF output from decompressing `test.cif.fcz` with format="mmcif".

---

### Directory and Archive Test Inputs

#### `test/dir_test_input/`
Directory containing copies of test structures for batch processing tests:
- `multichain.pdb` (4,678 atoms, chains A,B)
- `test.pdb` (2,208 atoms, chain A)
- `test.cif.gz` (1,513 atoms, chain A)
- `test_af.pdb` (243 atoms, chain A)

#### `test/tar_test_input.tar`
Uncompressed tar archive of test inputs.

#### `test/gz_test_input.tar.gz`
Gzipped tar archive of test inputs.
