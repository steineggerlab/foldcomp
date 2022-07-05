/**
 * File: main.cpp
 * Project: FoldU_background_search
 * Created: 2021-12-23 17:44:53
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     This code is written as part of project "src".
 * Usage:
 *    foldcomp compress input.pdb output.fcz
 *    foldcomp decompress input.fcz output.pdb
 * ---
 * Last Modified: 2022-06-27 11:11:14
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

int print_usage_with_error(void) {
    std::cout << "Usage: foldcomp compress input.pdb" << std::endl;
    std::cout << "       foldcomp compress input.pdb output.fcz" << std::endl;
    std::cout << "       foldcomp compress pdb_dir" << std::endl;
    std::cout << "       foldcomp compress pdb_dir fcz_dir" << std::endl;
    std::cout << "       foldcomp decompress input.fcz" << std::endl;
    std::cout << "       foldcomp decompress input.fcz output.pdb" << std::endl;
    std::cout << "       foldcomp decompress fcz_dir" << std::endl;
    std::cout << "       foldcomp decompress fcz_dir pdb_dir" << std::endl;
    return 1;
}

int compress(std::string input, std::string output) {
    // Check file extension
    if (input.substr(input.length() - 4) != ".pdb") {
        std::cerr << "Input file must be a .pdb file." << std::endl;
        return 1;
    }
    PDBReader pdbReader;
    pdbReader.load(input);
    std::vector<AtomCoordinate> atomCoordinates;
    atomCoordinates = pdbReader.loadAllAtomCoordinatesWithSideChain();
    std::string title = pdbReader.readTitle();
    std::vector<BackboneChain> compData;
    CompressedResidue compRes = CompressedResidue();
    // Convert title to char
    char* title_c = new char[title.length() + 1];
    strcpy(title_c, title.c_str());
    compRes.title = title_c;
    compRes.lenTitle = title.length() + 1;
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
    writeAtomCoordinatesToPDB(atomCoordinates, output);
    // DEBUGGING
    Nerf nerf;
    // nerf.writeInfoForChecking(atomCoordinates, "AFTER_COMPRESSION.csv");
    return 0;
}

int main(int argc, char const *argv[]) {
    if (argc < 3) {
        return print_usage_with_error();
    }
    enum {
        COMPRESS,
        DECOMPRESS,
        COMPRESS_MULTIPLE,
        DECOMPRESS_MULTIPLE
    } mode = COMPRESS;


    struct stat st = {0};
    if (stat(argv[2], &st) == -1) {
        std::cerr << "Error: " << argv[2] << " does not exist." << std::endl;
        return 1;
    }

    // get mode from command line
    if (strcmp(argv[1], "compress") == 0) {
        // Check argv[2] is file or directory
        // If directory, mode = COMPRESS_MULTIPLE
        // If file, mode = COMPRESS
        if (st.st_mode & S_ISDIR(st.st_mode)) {
            mode = COMPRESS_MULTIPLE;
        } else {
            mode = COMPRESS;
        }
    } else if (strcmp(argv[1], "decompress") == 0) {
        // Check argv[2] is file or directory
        // If directory, mode = DECOMPRESS_MULTIPLE
        // If file, mode = DECOMPRESS
        if (st.st_mode & S_ISDIR(st.st_mode)) {
            mode = DECOMPRESS_MULTIPLE;
        } else {
            mode = DECOMPRESS;
        }
    } else {
        return print_usage_with_error();
    }

    std::string input = argv[2];
    std::string output;
    if (argc == 4) {
        output = argv[3];
    }
    int flag;

    // check if mode is compress or decompress
    if (mode == COMPRESS) {
        // compress
        if (argc==3) {
            output = input.substr(0, input.length() - 4) + ".fcz";
        }
        std::cout << "Compressing " << input << " to " << output << std::endl;
        flag = compress(input, output);
    } else if (mode == DECOMPRESS) {
        // decompress
        if (argc==3) {
            output = input.substr(0, input.length() - 4) + ".pdb";
        }
        std::cout << "Decompressing " << input << " to " << output << std::endl;
        flag = decompress(input, output);
    } else if (mode == COMPRESS_MULTIPLE) {
        // compress multiple files
        // Check argument count
        if (argc == 3) {
            if (input[input.length() - 1] != '/') {
                input += "/";
            }
            output = input.substr(0, input.length() - 1) + "_fcz/";
        }
        if (output[output.length() - 1] != '/') {
            output += "/";
        }
        // Check output directory exists or not
        if (stat(output.c_str(), &st) == -1) {
        #if defined(_WIN32) || defined(_WIN64)
            flag = _mkdir(output.c_str());
        #else
            flag = mkdir(output.c_str(), 0777);
        #endif
        }
        // Get all files in input directory
        std::string file;
        std::string inputFile;
        std::string outputFile;
        std::vector<std::string> files = getFilesInDirectory(input);

        for (std::string file : files) {
            // Check file extension
            if (file.substr(file.length() - 4) == ".pdb") {
                outputFile = output + file.substr(0, file.length() - 4) + ".fcz";
                // Compress file
                inputFile = input + file;
                flag = compress(inputFile, outputFile);
                if (flag != 0) {
                    std::cerr << "Error compressing " << file << std::endl;
                    return 1;
                }
            }
        }
    } else if (mode == DECOMPRESS_MULTIPLE) {
        // decompress multiple files
        // Check argument count and input directory exists or not
        if (argc == 3) {
            if (input[input.length() - 1] != '/') {
                input += "/";
            }
            output = input.substr(0, input.length() - 1) + "_pdb/";
        }
        if (output[output.length() - 1] != '/') {
            output += "/";
        }
        // Check output directory exists or not
        if (stat(output.c_str(), &st) == -1) {
        #if defined(_WIN32) || defined(_WIN64)
            flag = _mkdir(output.c_str());
        #else
            flag = mkdir(output.c_str(), 0777);
        #endif
        }
        // Get all files in input directory
        std::string file;
        std::string inputFile;
        std::string outputFile;
        std::cout << "Decompressing files in " << input << std::endl;
        std::cout << "Output directory: " << output << std::endl;
        std::vector<std::string> files = getFilesInDirectory(input);
        for (std::string file : files) {
            // Check file extension
            if (file.substr(file.length() - 4) == ".fcz") {
                inputFile = input + file;
                outputFile = output + file.substr(0, file.length() - 4) + ".pdb";
                // Decompress file
                flag = decompress(inputFile, outputFile);
                if (flag != 0) {
                    std::cerr << "Error decompressing " << file << std::endl;
                    return 1;
                }
            }
        }
    } else {
        std::cout << "Invalid mode." << std::endl;
        return 1;
    }
    // Print log
    if (flag) {
        std::cout << "Error occurred." << std::endl;
    } else {
        std::cout << "Done." << std::endl;
    }
    return flag;
}
