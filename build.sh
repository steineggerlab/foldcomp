# File: build.sh
# Project: foldcomp
# Created: 2022-05-30 17:11:11
# Author: Hyunbin Kim (khb7840@gmail.com)
# Description:
#     Build script for foldcomp.
# ---
# Last Modified: 2022-09-20 11:52:22
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
./build/foldcomp compress test.pdb compressed.fcz
./build/foldcomp decompress compressed.fcz decompressed.pdb
# Test - Gzipped CIF
./build/foldcomp compress test.cif.gz compressed_cif.fcz
./build/foldcomp decompress compressed_cif.fcz decompressed_cif.pdb