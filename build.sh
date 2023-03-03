#!/bin/sh -e
# File: build.sh
# Project: foldcomp
# Created: 2022-05-30 17:11:11
# Author: Hyunbin Kim (khb7840@gmail.com)
# Description:
#     Build script for foldcomp.
# ---
# Last Modified: Fri Mar 03 2023
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
    ./build/foldcomp compress -y ./test/test.pdb ./test/compressed.fcz
    ./build/foldcomp decompress -y ./test/compressed.fcz ./test/decompressed.pdb
    # Test - Gzipped CIF
    ./build/foldcomp compress -y ./test/test.cif.gz ./test/compressed_cif.fcz
    ./build/foldcomp decompress -y -a ./test/compressed_cif.fcz ./test/decompressed_cif.pdb
    # RMSD
    RMSD1=$(./build/foldcomp rmsd ./test/test.pdb ./test/decompressed.pdb | cut -f6)
    awk -v check=$RMSD1 -v target=0.0826751 'BEGIN { diff = check - target; if (diff < 0) diff = -diff; if (diff > 0.001) { print check"!="target; exit 1 }  }'
    RMSD2=$(./build/foldcomp rmsd ./test/test.cif.gz ./test/decompressed_cif.pdb | cut -f6)
    awk -v check=$RMSD2 -v target=0.130284 'BEGIN { diff = check - target; if (diff < 0) diff = -diff; if (diff > 0.001) { print check"!="target; exit 1 }  }'
}

input_type_test()
{
    foldcomp_compress="./build/foldcomp compress -t 2 -y"
    foldcomp_decompress="./build/foldcomp decompress -t 2 -y"
    dir_input="./test/dir_test_input"
    tar_input="./test/tar_test_input.tar"
    gz_input="./test/gz_test_input.tar.gz"
    # TODO: implement database input with given ids
    # db_input="./test/db_test_input/pdb_db"
    # fcz_db_input="./test/example_db"

    # Compression - Directory
    echo "[Test1] Compression - Directory"
    # 01. Output is not given
    $foldcomp_compress $dir_input
    # 02. Output is given as
    # 02-1. Directory
    $foldcomp_compress $dir_input ./test/out/dir_in_dir_out
    # 02-2. Tarball
    $foldcomp_compress --tar $dir_input ./test/out/dir_in_tar_out.fcz.tar
    # 02-3. Database
    $foldcomp_compress --db $dir_input ./test/out/dir_in_db_out_fcz_db

    # Compression - Tarball
    echo "[Test2] Compression - Tarball"
    # 01. Output is not given
    $foldcomp_compress $tar_input
    # 02. Output is given as
    # 02-1. Directory
    $foldcomp_compress $tar_input ./test/out/tar_in_dir_out
    # 02-2. Tarball
    $foldcomp_compress --tar $tar_input ./test/out/tar_in_tar_out.fcz.tar
    # 02-3. Database
    $foldcomp_compress --db $tar_input ./test/out/tar_in_db_out_fcz_db

    # Compression - Gzipped Tarball
    echo "[Test3] Compression - Gzipped Tarball"
    # 01. Output is not given
    $foldcomp_compress $gz_input
    # 02. Output is given as
    # 02-1. Directory
    $foldcomp_compress $gz_input ./test/out/gz_in_dir_out
    # 02-2. Tarball
    $foldcomp_compress --tar $gz_input ./test/out/gz_in_tar_out.fcz.tar
    # 02-3. Database
    $foldcomp_compress --db $gz_input ./test/out/gz_in_db_out_fcz_db

    # # Compression - Database
    # echo "[Test4] Compression - Database"
    # # 01. Output is not given
    # $foldcomp_compress $db_input
    # # 02. Output is given as
    # # 02-1. Directory
    # $foldcomp_compress $db_input ./test/out/db_in_dir_out
    # # 02-2. Tarball
    # $foldcomp_compress --tar $db_input ./test/out/db_in_tar_out.tar
    # # 02-3. Database
    # $foldcomp_compress --db $db_input ./test/out/db_in_db_out

    # Decompression - directory
    echo "[Test5] Decompression - Directory"
    # 01. Output is not given
    $foldcomp_decompress ${dir_input}_fcz
    # 02. Output is given as
    # 02-1. Directory
    $foldcomp_decompress ${dir_input}_fcz ./test/out/dir_in_dir_out
    # 02-2. Tarball
    $foldcomp_decompress --tar ${dir_input}_fcz ./test/out/dir_in_tar_out.pdb.tar
    # 02-3. Database
    $foldcomp_decompress --db ${dir_input}_fcz ./test/out/dir_in_db_out_pdb_db

    # Decompression - Tarball
    echo "[Test6] Decompression - Tarball"
    # 01. Output is not given
    $foldcomp_decompress ./test/out/dir_in_tar_out.fcz.tar
    # 02. Output is given as
    # 02-1. Directory
    $foldcomp_decompress ./test/out/dir_in_tar_out.fcz.tar ./test/out/tar_in_dir_out
    # 02-2. Tarball
    $foldcomp_decompress --tar ./test/out/dir_in_tar_out.fcz.tar ./test/out/tar_in_tar_out.pdb.tar
    # 02-3. Database
    $foldcomp_decompress --db ./test/out/dir_in_tar_out.fcz.tar ./test/out/tar_in_db_out_pdb_db

    # Decompression - Gzipped Tarball
    echo "[Test7] Decompression - Gzipped Tarball"
    # 01. Output is not given
    gzip -f ./test/out/dir_in_tar_out.fcz.tar
    $foldcomp_decompress ./test/out/dir_in_tar_out.fcz.tar.gz
    # 02. Output is given as
    # 02-1. Directory
    $foldcomp_decompress ./test/out/dir_in_tar_out.fcz.tar.gz ./test/out/gz_in_dir_out
    # 02-2. Tarball
    $foldcomp_decompress --tar ./test/out/dir_in_tar_out.fcz.tar.gz ./test/out/gz_in_tar_out.pdb.tar
    # 02-3. Database
    $foldcomp_decompress --db ./test/out/dir_in_tar_out.fcz.tar.gz ./test/out/gz_in_db_out_pdb_db

    # # Decompression - Database
    # echo "[Test8] Decompression - Database"
    # # 01. Output is not given
    # $foldcomp_decompress $fcz_db_input
    # # 02. Output is given as
    # # 02-1. Directory
    # $foldcomp_decompress $fcz_db_input ./test/out/db_in_dir_out
    # # 02-2. Tarball
    # $foldcomp_decompress --tar $fcz_db_input ./test/out/db_in_tar_out.pdb.tar
    # # 02-3. Database
    # $foldcomp_decompress --db $fcz_db_input ./test/out/db_in_db_out_pdb_db
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