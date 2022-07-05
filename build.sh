# File: build.sh
# Project: FoldU_background_search
# Created: 2022-05-30 17:11:11
# Author: Hyunbin Kim (khb7840@gmail.com)
# Description:
#     This code is written as part of project "FoldU_background_search".
# ---
# Last Modified: 2022-05-31 18:05:46
# Modified By: Hyunbin Kim (khb7840@gmail.com)
# ---
# Copyright © 2022 Hyunbin Kim, All rights reserved

# Configure
mkdir build
cd build
cmake ..
cd ..
# Build
cmake --build ./build --target foldcomp

# Test
./build/foldcomp compress test.pdb compressed.fcz
./build/foldcomp decompress compressed.fcz decompressed.pdb