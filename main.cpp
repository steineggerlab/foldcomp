/**
 * File: main.cpp
 * Project: hbk
 * Created: 2021-12-23 17:44:53
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     This code contains main function for "foldcomp".
 *     Foldcomp is a fast lossy compression algorithm for protein structure.
 *     It encodes torsion angles with optimal number of bits and reconstruct
 *     3D coordinates from the encoded angles.
 * Usage:
 *    foldcomp compress input.pdb output.fcz
 *    foldcomp decompress input.fcz output.pdb
 * ---
 * Last Modified: 2022-07-20 06:53:58
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
// Headers in the project
#include "aa_sidechain.h"
#include "amino_acid.h"
#include "atom_coordinate.h"
#include "compressed_torsion.h"
#include "discretizer.h"
#include "nerf.h"
#include "pdb_reader.h"
#include "torsion.xyz.h"
#include "utility.h"
// Standard libraries
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <sys/stat.h>
#include <cstring>
#include <getopt.h>
// OpenMP for parallelization
#include <omp.h>

int print_usage(void) {
    std::cout << "Usage: foldcomp compress <pdb_file> [<fcz_file>]" << std::endl;
    std::cout << "       foldcomp compress [-t number] <pdb_dir> [<fcz_dir>]" << std::endl;
    std::cout << "       foldcomp decompress <fcz_file> [<pdb_file>]" << std::endl;
    std::cout << "       foldcomp decompress [-t number] <fcz_dir> [<pdb_dir>]" << std::endl;
    std::cout << " -t, --threads        number of threads to use [default=1]" << std::endl;
    std::cout << " -h, --help           print this help message" << std::endl;
    return 1;
}

int compress(std::string input, std::string output) {
    PDBReader pdbReader;
    pdbReader.load(input);
    std::vector<AtomCoordinate> atomCoordinates;
    atomCoordinates = pdbReader.loadAllAtomCoordinatesWithSideChain();
    std::string title = pdbReader.readTitle();
    pdbReader.infile.close();
    std::vector<BackboneChain> compData;
    CompressedResidue compRes = CompressedResidue();
    // Convert title to char
    compRes.strTitle = title;
    compData = compRes.compress(atomCoordinates);
    // Write compressed data to file
    compRes.write(output);
    // DEBUGGING
    Nerf nerf;
    // nerf.writeInfoForChecking(atomCoordinates, "BEFORE_COMPRESSION.csv");
    // clear memory
    atomCoordinates.clear();
    compData.clear();
    return 0;
}

int decompress(std::string input, std::string output) {
    int flag = 0;
    CompressedResidue compRes = CompressedResidue();
    flag = compRes.read(input);
    if (flag != 0) {
        std::cerr << "Error reading " << input << std::endl;
        return 1;
    }
    std::vector<AtomCoordinate> atomCoordinates;
    flag = compRes.decompress(atomCoordinates);
    if (flag != 0) {
        std::cerr << "Error decompressing compressed data." << std::endl;
        return 1;
    }
    // Write decompressed data to file
    flag = writeAtomCoordinatesToPDB(atomCoordinates, compRes.strTitle, output);

    return flag;
}

std::string getFileWithoutExt(std::string& file) {
    size_t extStart = file.find_last_of('.');
    return extStart == std::string::npos ? file : file.substr(0, extStart);
}

int main(int argc, char* const *argv) {
    if (argc < 3) {
        return print_usage();
    }
    
    int flag = 0;
    int option_index = 0;
    int num_threads = 1;
    int has_output = 0;

    // Mode - non-optional argument
    enum {
        COMPRESS,
        DECOMPRESS,
        COMPRESS_MULTIPLE,
        DECOMPRESS_MULTIPLE
    } mode = COMPRESS;

    // Define command line options
    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"threads", required_argument, 0, 't'},
        {0, 0, 0, 0}
    };

    // Parse command line options with getopt_long
    flag = getopt_long(argc, argv, "ht:", long_options, &option_index);

    while (flag != -1) {
        switch (flag) {
        case 'h':
            return print_usage();
        case 't':
            num_threads = atoi(optarg);
            break;
        case '?':
            return print_usage();
        default:
            break;
        }
        flag = getopt_long(argc, argv, "ht:", long_options, &option_index);
    }
    // // Set number of threads
    // omp_set_num_threads(num_threads);

    // Parse non-option arguments
    // argv[optind]: MODE
    // argv[optind + 1]: INPUT
    // argv[optind + 2]: OUTPUT (optional)

    if ((optind + 1) >= argc) {
        std::cerr << "Error: Not enough arguments." << std::endl;
        return print_usage();
    }

    struct stat st = {0};
    if (stat(argv[optind + 1], &st) == -1) {
        std::cerr << "Error: " << argv[optind + 1] << " does not exist." << std::endl;
        return 1;
    }

    // get mode from command line
    if (strcmp(argv[optind], "compress") == 0) {
        // Check argv[2] is file or directory
        // If directory, mode = COMPRESS_MULTIPLE
        // If file, mode = COMPRESS
        if (st.st_mode & S_ISDIR(st.st_mode)) {
            mode = COMPRESS_MULTIPLE;
        } else {
            mode = COMPRESS;
        }
    } else if (strcmp(argv[optind], "decompress") == 0) {
        // Check argv[2] is file or directory
        // If directory, mode = DECOMPRESS_MULTIPLE
        // If file, mode = DECOMPRESS
        if (st.st_mode & S_ISDIR(st.st_mode)) {
            mode = DECOMPRESS_MULTIPLE;
        } else {
            mode = DECOMPRESS;
        }
    } else {
        return print_usage();
    }
    std::string input = argv[optind + 1];
    std::string output;
    if (argc == optind + 3) {
        has_output = 1;
        output = argv[optind + 2];
    }

    // check if mode is compress or decompress
    if (mode == COMPRESS) {
        // compress
        if (!has_output) {
            output = getFileWithoutExt(input) + ".fcz";
        }
        std::cout << "Compressing " << input << " to " << output << std::endl;
        compress(input, output);
        flag = 0;
    } else if (mode == DECOMPRESS) {
        // decompress
        if (!has_output) {
            output = getFileWithoutExt(input) + "_fcz.pdb";
        }
        std::cout << "Decompressing " << input << " to " << output << std::endl;
        decompress(input, output);
        flag = 0;
    } else if (mode == COMPRESS_MULTIPLE) {
        // compress multiple files
        if (input[input.length() - 1] != '/') {
            input += "/";
        }
        if (!has_output) {
            output = input.substr(0, input.length() - 1) + "_fcz/";
        }
        if (output[output.length() - 1] != '/') {
            output += "/";
        }
        // Check output directory exists or not
        if (stat(output.c_str(), &st) == -1) {
        #if defined(_WIN32) || defined(_WIN64)
            _mkdir(output.c_str());
        #else
            mkdir(output.c_str(), 0777);
        #endif
        }
        // Get all files in input directory
        std::string file;
        std::string inputFile;
        std::string outputFile;
        std::cout << "Compressing files in " << input;
        std::cout << " using " << num_threads << " threads" <<std::endl;
        std::cout << "Output directory: " << output << std::endl;
        std::vector<std::string> files = getFilesInDirectory(input);
        // Parallelize
        omp_set_num_threads(num_threads);
        #pragma omp parallel
        {
            #pragma omp for
            for (int i = 0; i < files.size(); i++) {
                std::string file = files[i];
                std::string inputFile = input + file;
                std::string outputFile = output + getFileWithoutExt(file) + ".fcz";
                compress(inputFile, outputFile);
            }
        }
        flag = 0;
    } else if (mode == DECOMPRESS_MULTIPLE) {
        // decompress multiple files
        if (input[input.length() - 1] != '/') {
            input += "/";
        }
        if (!has_output) {
            output = input.substr(0, input.length() - 1) + "_pdb/";
        }
        if (output[output.length() - 1] != '/') {
            output += "/";
        }
        // Check output directory exists or not
        if (stat(output.c_str(), &st) == -1) {
        #if defined(_WIN32) || defined(_WIN64)
            _mkdir(output.c_str());
        #else
            mkdir(output.c_str(), 0777);
        #endif
        }
        // Get all files in input directory
        std::string file;
        std::string inputFile;
        std::string outputFile;
        std::cout << "Decompressing files in " << input;
        std::cout << " using " << num_threads << " threads" << std::endl;
        std::cout << "Output directory: " << output << std::endl;
        std::vector<std::string> files = getFilesInDirectory(input);
        omp_set_num_threads(num_threads);
        #pragma omp parallel
        {
            #pragma omp for
            for (int i = 0; i < files.size(); i++) {
                std::string file = files[i];
                std::string inputFile = input + file;
                std::string outputFile = output + file.substr(0, file.length() - 4) + ".pdb";
                decompress(inputFile, outputFile);
            }
        }
        flag = 0;
    } else {
        std::cout << "Invalid mode." << std::endl;
        return 1;
    }
    // Print log
    std::cout << "Done." << std::endl;
    return flag;
}
