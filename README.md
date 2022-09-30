# Foldcomp
Foldcomp compresses protein structures with torsion angles effectively. It compresses the backbone atoms to 8 bytes and the side chain to additionally 4-5 byes per residue, an averaged-sized protein of 350 residues requires ~4.2kb.

## Usage

### Executable
```
[Compression]
foldcomp compress <pdb_file|cif_file> [<fcz_file>]
foldcomp compress [-t number] <pdb_dir|cif_dir> [<fcz_dir>]

[Decompression]
foldcomp decompress <fcz_file> [<pdb_file>]
foldcomp decompress [-t number] <fcz_dir> [<pdb_dir>]

[Extraction]
foldcomp extract [--plddt|--fasta] <fcz_file> [<txt_file|fasta_file>]
foldcomp extract [--plddt|--fasta] [-t number] <fcz_dir|tar> [<output_dir>]

[Check]
foldcomp check <fcz_file>
foldcomp check [-t number] <fcz_dir|tar>

[Options]
 -h, --help           print this help message
 -t, --threads        number of threads to use [default=1]
 -a, --alt            use alternative atom order [default=false]
 -b, --break          interval size to save absolute atom coordinates [default=200]
 -z, --tar            save as tar file [default=false]
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
fcz = open("test/compressed.fcz", "rb")
fcz_binary = fcz.read()
fcz.close()

# Decompress
pdb_out = foldcomp.decompress(fcz_binary) # pdb_out[0]: file name, pdb_out[1]: pdb binary string

# Save to a pdb file
out_file = open(pdb_out[0], "wb")
out_file.write(pdb_out[1])
out_file.close()

# 02. Iterate over a database of FCZ files
# Open a foldcomp database
example_db = foldcomp.open(path=Path("test/example_db"))
ids = ["d1asha_", "d1it2a_"]

# Iterate through database
for entry in example_db:
    # entry[0]: id, entry[1]: pdb binary string
    if entry[0] in ids:
        # Save matched entries to pdb files
        out_file = open(entry[0] + ".pdb", "wb")
        out_file.write(entry[1])
        out_file.close()
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

Foldcomp is a compression method and format to compress protein structures requiring only 13 bytes per residue, which reduces the required storage space by an order of magnitude than saving 3D coordinates directly. We achieve this reduction by encoding the torsion angles of the backbone as well as the side-chain angles in a compact binary file format, FCZ.
We

![abstract](.github/img/Abstract.jpg)

> WARNING: Current version of Foldcomp does not support compression of multiple chains in a single file.

## Contributor
<a href="https://github.com/steineggerlab/foldcomp/graphs/contributors">
  <img src="https://contributors-img.firebaseapp.com/image?repo=steineggerlab/foldcomp" />
</a>
