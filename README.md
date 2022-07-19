---
File: README.md
ProjectName: Foldcomp
Author: Hyunbin Kim (khb7840@gmail.com)
Last Modified: 2022-07-20 02:02:07
---

# Foldcomp
Foldcomp compresses protein structures with torsion angles effectively.

## Usage

### Build
```sh
# Configure
mkdir build
cd build
cmake ..
cd ..
# Build
cmake --build ./build --target foldcomp
```
### Compression
```sh
# Single PDB file
./build/foldcomp compress test.pdb  # default output: test.fcz
./build/foldcomp compress test.pdb compressed.fcz

# Directory of PDB files
./build/foldcomp compress pdb_dir/ # default output: pdb_dir_fcz/
./build/foldcomp compress pdb_dir/ compressed_dir/
```
### Decompression
```sh
# Single FCZ file
./build/foldcomp decompress compressed.fcz # default output: compressed.pdb
./build/foldcomp decompress compressed.fcz decompressed.pdb

# Directory of FCZ files
./build/foldcomp decompress fcz_dir/ # default output: fcz_dir_pdb/
./build/foldcomp decompress fcz_dir/ decompressed_dir/
```

---

![abstract](.github/img/Abstract.jpg)