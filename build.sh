#!/bin/sh -e
# File: build.sh
# Project: foldcomp
# Created: 2022-05-30 17:11:11
# Author: Hyunbin Kim (khb7840@gmail.com)
# Description:
#     Build script for foldcomp.
# ---
# Last Modified: Tue Feb 28 2023
# Modified By: Hyunbin Kim
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

minimal_test()
{
    # Test - PDB
    ./build/foldcomp compress ./test/test.pdb ./test/compressed.fcz
    ./build/foldcomp decompress ./test/compressed.fcz ./test/decompressed.pdb
    # Test - Gzipped CIF
    ./build/foldcomp compress ./test/test.cif.gz ./test/compressed_cif.fcz
    ./build/foldcomp decompress -a ./test/compressed_cif.fcz ./test/decompressed_cif.pdb
    # RMSD
    RMSD1=$(./build/foldcomp rmsd ./test/test.pdb ./test/decompressed.pdb | cut -f6)
    awk -v check=$RMSD1 -v target=0.0826751 'BEGIN { diff = check - target; if (diff < 0) diff = -diff; if (diff > 0.001) { print check"!="target; exit 1 }  }'
    RMSD2=$(./build/foldcomp rmsd ./test/test.cif.gz ./test/decompressed_cif.pdb | cut -f6)
    awk -v check=$RMSD2 -v target=0.130284 'BEGIN { diff = check - target; if (diff < 0) diff = -diff; if (diff > 0.001) { print check"!="target; exit 1 }  }'
}

input_type_test()
{
    dir_input="./test/dir_test_input"
    tar_input="./test/tar_test_input.tar"
    gz_input="./test/gz_test_input.tar.gz"

    # Compression - Directory
    # 01. Output is not given
    ./build/foldcomp compress $dir_input
    # 02. Output is given as
    # 02-1. Directory
    ./build/foldcomp compress $dir_input ./test/out/dir_in_dir_out
    # 02-2. Tarball
    ./build/foldcomp compress --tar $dir_input ./test/out/dir_in_tar_out.tar
    # 02-3. Database
    ./build/foldcomp compress --db $dir_input ./test/out/dir_in_db_out

    # Compression - Tarball
    # 01. Output is not given
    ./build/foldcomp compress $tar_input
    # 02. Output is given as
    # 02-1. Directory
    ./build/foldcomp compress $tar_input ./test/out/tar_in_dir_out
    # 02-2. Tarball
    ./build/foldcomp compress --tar $tar_input ./test/out/tar_in_tar_out.tar
    # 02-3. Database
    ./build/foldcomp compress --db $tar_input ./test/out/tar_in_db_out

    # Compression - Gzipped Tarball
    # 01. Output is not given
    ./build/foldcomp compress $gz_input
    # 02. Output is given as
    # 02-1. Directory
    ./build/foldcomp compress $gz_input ./test/out/gz_in_dir_out
    # 02-2. Tarball
    ./build/foldcomp compress --tar $gz_input ./test/out/gz_in_tar_out.tar
    # 02-3. Database
    ./build/foldcomp compress --db $gz_input ./test/out/gz_in_db_out

    # Decompression - database
    # TODO: Implement
}


# Check argument is not given
if [ $# -eq 0 ]; then
    # Run minimal test
    minimal_test
    echo "All good!"
    exit 0
else
    # Check argument is "test"
    if [ $1 = "test" ]; then
        # Run minimal test
        minimal_test
        # Run full test
        input_type_test
        echo "All good!"
        exit 0
    else
        echo "Invalid argument"
        exit 1
    fi
fi