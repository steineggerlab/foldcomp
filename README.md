# Foldcomp
Foldcomp compresses protein structures with torsion angles effectively. It compresses the backbone atoms to 8 bytes and the side chain to additionally 4-5 byes per residue, thus an averaged-sized protein of 350 residues requires ~4.2kb.

## Usage

### Executable
```
# Compression
foldcomp compress <pdb_file|cif_file> [<fcz_file>]
foldcomp compress [-t number] <pdb_dir|cif_dir> [<fcz_dir>]

# Decompression
foldcomp decompress <fcz_file> [<pdb_file>]
foldcomp decompress [-t number] <fcz_dir> [<pdb_dir>]

# Extraction of sequence or pLDDT
foldcomp extract [--plddt|--fasta] <fcz_file> [<txt_file|fasta_file>]
foldcomp extract [--plddt|--fasta] [-t number] <fcz_dir|tar> [<output_dir>]

# Check
foldcomp check <fcz_file>
foldcomp check [-t number] <fcz_dir|tar>

# Options
 -h, --help           print this help message
 -t, --threads        number of threads to use [default=1]
 -a, --alt            use alternative atom order [default=false]
 -b, --break          interval size to save absolute atom coordinates [default=200]
 -z, --tar            save as tar file [default=false]
 -c, --cpu            CPU cores for (de)compression of folders/tar files [default=1]
 --plddt              extract pLDDT score (only for extraction mode)
 --fasta              extract amino acid sequence (only for extraction mode)
 --no-merge           do not merge output files (only for extraction mode)
```

### Python API
```py
import foldcomp
from pathlib import Path
# 01. Handling a FCZ file
# Open a fcz file
with open("test/compressed.fcz", "rb") as fcz:
  fcz_binary = fcz.read()

  # Decompress
  (name, pdb) = foldcomp.decompress(fcz_binary) # pdb_out[0]: file name, pdb_out[1]: pdb binary string

  # Save to a pdb file
  with open(name, "wb") as pdb_file:
    pdb_file.write(pdb)

# 02. Iterate over a database of FCZ files
# Open a foldcomp database
ids = ["d1asha_", "d1it2a_"]
with foldcomp.open(path=Path("test/example_db"), uniprot_ids=ids) as db:
  # Iterate through database
  for (name, pdb) in db:
      # save entries as seperate pdb files
      with open(name + ".pdb", "wb") as pdb_file:
        pdb_file.write(pdb)
```

## Build

### Executable
```sh
./build.sh
```
or
```sh
# Configure
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cd ..
# Build
cmake --build ./build --target foldcomp
```

## About

Foldcomp is a compression method and format to compress protein structures requiring only 13 bytes per residue, which reduces the required storage space by an order of magnitude compared to saving 3D coordinates directly. We achieve this reduction by encoding the torsion angles of the backbone as well as the side-chain angles in a compact binary file format, FCZ.

![abstract](.github/img/Abstract.jpg)

> WARNING: Current version of Foldcomp does not support compression of multiple chains in a single file.

## Contributor
<a href="https://github.com/steineggerlab/foldcomp/graphs/contributors">
  <img src="https://contributors-img.firebaseapp.com/image?repo=steineggerlab/foldcomp" />
</a>
