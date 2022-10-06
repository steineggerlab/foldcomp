# File: build.sh
# Project: foldcomp
# Created: 2022-05-30 17:11:11
# Author: Hyunbin Kim (khb7840@gmail.com)
# Description:
#     Build script for foldcomp.
# ---
# Last Modified: 2022-10-06 20:38:41
# Modified By: Hyunbin Kim (khb7840@gmail.com)
# ---
# Copyright © 2022 Hyunbin Kim, All rights reserved

# Configure
# If directory "build" does not exist, create it.
if [ ! -d "build" ]; then
    mkdir build
fi
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cd ..
# Build
cmake --build ./build --target foldcomp

# Test - PDB
./build/foldcomp compress ./test/test.pdb ./test/compressed.fcz
./build/foldcomp decompress ./test/compressed.fcz ./test/decompressed.pdb
./build/foldcomp rmsd ./test/test.pdb ./test/decompressed.pdb
# Test - Gzipped CIF
./build/foldcomp compress ./test/test.cif.gz ./test/compressed_cif.fcz
./build/foldcomp decompress -a ./test/compressed_cif.fcz ./test/decompressed_cif.pdb