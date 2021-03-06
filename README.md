---
File: README.md
ProjectName: Foldcomp
Author: Hyunbin Kim (khb7840@gmail.com)
Last Modified: 2022-07-20 12:31:51
---

# Foldcomp
Foldcomp compresses protein structures with torsion angles effectively.

## Usage
```
[Compression]
foldcomp compress <pdb_file> [<fcz_file>]
foldcomp compress [-t number] <pdb_dir> [<fcz_dir>]

[Decompression]
foldcomp decompress <fcz_file> [<pdb_file>]
foldcomp decompress [-t number] <fcz_dir> [<pdb_dir>]

[Options]
 -h, --help           print this help message
 -t, --threads        number of threads to use [default=1]
 -a, --alt            use alternative atom order [default=false]
```

## Build
```sh
# Configure
mkdir build
cd build
cmake ..
cd ..
# Build
cmake --build ./build --target foldcomp
```


## About

Foldcomp is a compression method and format to compress protein structures requiring only 13 bytes per residue, which reduces the required storage space by an order of magnitude than saving 3D coordinates directly. We achieve this reduction by encoding the torsion angles of the backbone as well as the side-chain angles in a compact format.

![abstract](.github/img/Abstract.jpg)
