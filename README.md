# Foldcomp
Foldcomp compresses protein structures with torsion angles effectively. It compresses the backbone atoms to 8 bytes and the side chain to additionally 4-5 byes pre residue, an averaged sized protein of 350 residues require ~4.2kb.


## Usage
```
[Compression]
foldcomp compress <pdb_file|cif_file> [<fcz_file>]
foldcomp compress [-t number] <pdb_dir|cif_dir> [<fcz_dir>]

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
cmake -DCMAKE_BUILD_TYPE=Release ..
cd ..
# Build
cmake --build ./build --target foldcomp
```


## About

Foldcomp is a compression method and format to compress protein structures requiring only 13 bytes per residue, which reduces the required storage space by an order of magnitude than saving 3D coordinates directly. We achieve this reduction by encoding the torsion angles of the backbone as well as the side-chain angles in a compact format.

![abstract](.github/img/Abstract.jpg)

> WARNING: Current version of Foldcomp does not support compression of multiple chains in a single file.

## Contributor
<a href="https://github.com/steineggerlab/foldcomp/graphs/contributors">
  <img src="https://contributors-img.firebaseapp.com/image?repo=steineggerlab/foldcomp" />
</a>
