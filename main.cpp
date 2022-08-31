/**
 * File: main.cpp
 * Project: foldcomp
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
 * Last Modified: 2022-08-31 13:34:46
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
#include "input.h"
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

#ifdef HAVE_GCS
#include "google/cloud/storage/client.h"
#endif

static int use_alt_order = 0;
static int anchor_residue_threshold = 200;
static int save_as_tar = 0;

int print_usage(void) {
    std::cout << "Usage: foldcomp compress <pdb_file> [<fcz_file>]" << std::endl;
    std::cout << "       foldcomp compress [-t number] <pdb_dir> [<fcz_dir>]" << std::endl;
    std::cout << "       foldcomp decompress <fcz_file> [<pdb_file>]" << std::endl;
    std::cout << "       foldcomp decompress [-t number] <fcz_dir> [<pdb_dir>]" << std::endl;
    std::cout << " -h, --help           print this help message" << std::endl;
    std::cout << " -t, --threads        number of threads to use [default=1]" << std::endl;
    std::cout << " -a, --alt            use alternative atom order [default=false]" << std::endl;
    std::cout << " -b, --break          interval size to save absolute atom coordinates [default=200]" << std::endl;
    std::cout << " -z, --tar                save as tar file [default=false]" << std::endl;
    return 0;
}

int compress(std::string input, std::string output) {
    StructureReader reader;
    reader.load(input);
    std::vector<AtomCoordinate> atomCoordinates;
    reader.readAllAtoms(atomCoordinates);
    if (atomCoordinates.size() == 0) {
        std::cout << "Error: No atoms found in the input file: " << input << std::endl;
        return 1;
    }
    std::string title = reader.title;

    std::vector<BackboneChain> compData;
    CompressedResidue compRes = CompressedResidue();
    // Convert title to char
    compRes.strTitle = title;
    compRes.anchorThreshold = anchor_residue_threshold;
    compData = compRes.compress(atomCoordinates);
    // Write compressed data to file
    compRes.write(output);
    // DEBUGGING
    // Nerf nerf;
    // nerf.writeInfoForChecking(atomCoordinates, "BEFORE_COMPRESSION.csv");
    // clear memory
    // atomCoordinates.clear();
    // compData.clear();
    return 0;
}


int compressFromBuffer(const std::string& content, const std::string& output, std::string& name) {
    StructureReader reader;
    reader.loadFromBuffer(content.c_str(), content.size(), name);
    std::vector<AtomCoordinate> atomCoordinates;
    reader.readAllAtoms(atomCoordinates);
    if (atomCoordinates.size() == 0) {
        std::cout << "Error: No atoms found in the input" << std::endl;
        return 1;
    }
    std::string title = reader.title;

    std::vector<BackboneChain> compData;
    CompressedResidue compRes = CompressedResidue();
    // Convert title to char
    compRes.strTitle = name;
    compRes.anchorThreshold = anchor_residue_threshold;
    compData = compRes.compress(atomCoordinates);
    // Write compressed data to file
    compRes.write(output + "/" + name + ".fcz");
    // DEBUGGING
    // Nerf nerf;
    // nerf.writeInfoForChecking(atomCoordinates, "BEFORE_COMPRESSION.csv");
    // clear memory
    // atomCoordinates.clear();
    // compData.clear();
    return 0;
}


int compressWithoutWriting(CompressedResidue& compRes, std::string input) {
    StructureReader reader;
    reader.load(input);
    std::vector<AtomCoordinate> atomCoordinates;
    reader.readAllAtoms(atomCoordinates);
    if (atomCoordinates.size() == 0) {
        std::cout << "Error: No atoms found in the input file: " << input << std::endl;
        return 1;
    }
    std::string title = reader.title;

    std::vector<BackboneChain> compData;
    // Convert title to char
    compRes.strTitle = title;
    compRes.anchorThreshold = anchor_residue_threshold;
    compData = compRes.compress(atomCoordinates);
    return 0;
}

int compressFromBufferWithoutWriting(CompressedResidue& compRes, const std::string& content, std::string& name) {
    StructureReader reader;
    reader.loadFromBuffer(content.c_str(), content.size(), name);
    std::vector<AtomCoordinate> atomCoordinates;
    reader.readAllAtoms(atomCoordinates);
    if (atomCoordinates.size() == 0) {
        std::cout << "Error: No atoms found in the input" << std::endl;
        return 1;
    }
    std::string title = name;
    std::vector<BackboneChain> compData;
    // Convert title to char
    compRes.strTitle = name;
    compRes.anchorThreshold = anchor_residue_threshold;
    compData = compRes.compress(atomCoordinates);
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
    compRes.useAltAtomOrder = use_alt_order;
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

bool stringEndsWith(const std::string& suffix, const std::string& str) {
    if (str.length() < suffix.length()) {
        return false;
    }

    return (!str.compare(str.length() - suffix.length(), suffix.length(), suffix));
}

bool stringStartsWith(const std::string& prefix, const std::string& str, const size_t offset = 0) {
    if (str.length() < prefix.length()) {
        return false;
    }
    return (!str.compare(offset, prefix.length(), prefix));
}


std::vector<std::string> stringSplit(const std::string& str, const std::string& sep) {
    std::vector<std::string> arr;

    char* cstr = strdup(str.c_str());
    const char* csep = sep.c_str();
    char* rest;
    char* current = strtok_r(cstr, csep, &rest);
    while (current != NULL) {
        arr.emplace_back(current);
        current = strtok_r(NULL, csep, &rest);
    }
    free(cstr);

    return arr;
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
        COMPRESS_MULTIPLE_GCS,
        DECOMPRESS_MULTIPLE
    } mode = COMPRESS;

    // Define command line options
    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"alt", no_argument, 0, 'a'},
        {"tar", no_argument, 0, 'z'},
        {"threads", required_argument, 0, 't'},
        {"break", required_argument, 0, 'b'},
        {0, 0, 0, 0}
    };

    // Parse command line options with getopt_long
    flag = getopt_long(argc, argv, "hazt:b:", long_options, &option_index);

    while (flag != -1) {
        switch (flag) {
        case 'h':
            return print_usage();
        case 't':
            num_threads = atoi(optarg);
            break;
        case 'a':
            use_alt_order = 1;
            break;
        case 'z':
            save_as_tar = 1;
            break;
        case 'b':
            anchor_residue_threshold = atoi(optarg);
            break;
        case '?':
            return print_usage();
        default:
            break;
        }
        flag = getopt_long(argc, argv, "hazt:b:", long_options, &option_index);
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

    struct stat st = { 0 };
    int fileExists = stat(argv[optind + 1], &st);
    // get mode from command line
    if (strcmp(argv[optind], "compress") == 0) {
        // Check argv[2] is file, directory, or gcs URI
        // TODO: stdin support
        // If gcs URI, mode = COMPRESS_MULTIPLE_GCS
        // If directory, mode = COMPRESS_MULTIPLE
        // If file, mode = COMPRESS
#ifdef HAVE_GCS
        if ((optind + 1) < argc && stringStartsWith("gcs://", argv[optind + 1])) {
            mode = COMPRESS_MULTIPLE_GCS;
            fileExists = 0;
        } else
#endif
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
    // Error if no input file given
    if (mode != COMPRESS_MULTIPLE_GCS && fileExists == -1) {
        std::cerr << "Error: " << argv[optind + 1] << " does not exist." << std::endl;
        return 1;
    }

    std::string input = argv[optind + 1];
    std::string output;
    if (argc == optind + 3) {
        has_output = 1;
        output = argv[optind + 2];
    }

    // check if mode is compress or decompress
    if (mode == COMPRESS) {
        // compress a single file
        if (!has_output) {
            output = getFileWithoutExt(input) + ".fcz";
        }
        std::cout << "Compressing " << input << " to " << output << std::endl;
        compress(input, output);
        flag = 0;
    } else if (mode == DECOMPRESS) {
        // decompress a single file
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
        } else {
            if (stringEndsWith(output, ".tar")) {
                save_as_tar = 1;
            }
        }

        if (output[output.length() - 1] != '/' && !save_as_tar) {
            output += "/";
        }
        // Check output directory exists or not
        if (!save_as_tar) {
            if (stat(output.c_str(), &st) == -1) {
            #if defined(_WIN32) || defined(_WIN64)
                _mkdir(output.c_str());
            #else
                mkdir(output.c_str(), 0755);
            #endif
            }
        }
        // Get all files in input directory
        std::string file;
        std::string inputFile;
        std::string outputFile;
        std::cout << "Compressing files in " << input;
        std::cout << " using " << num_threads << " threads" <<std::endl;
        if (save_as_tar) {
            std::cout << "Output tar file: " << output << std::endl;
        } else {
            std::cout << "Output directory: " << output << std::endl;
        }
        std::vector<std::string> files = getFilesInDirectory(input);
        // Parallelize
        omp_set_num_threads(num_threads);
        if (!save_as_tar) {
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
        } else {
            mtar_t tar;
            std::string tarFile = output.substr(0, output.length() - 1) + ".tar";
            mtar_open(&tar, tarFile.c_str(), "w");
            #pragma omp parallel
            {
            #pragma omp for
                for (int i = 0; i < files.size(); i++) {
                    std::string file = files[i];
                    std::string inputFile = input + file;
                    std::string outputFile = output + getFileWithoutExt(file) + ".fcz";
                    CompressedResidue compRes = CompressedResidue();
                    compressWithoutWriting(compRes, inputFile);
                    compRes.writeTar(tar, outputFile, compRes.getSize());
                }
            }
            mtar_finalize(&tar);
            mtar_close(&tar);
        }
    } else if (mode == COMPRESS_MULTIPLE_GCS) {
        // compress multiple files from gcs
#ifdef HAVE_GCS
        if (!has_output) {
            std::cerr << "Please specify output directory" << std::endl;
            return 1;
        }
        if (output[output.length() - 1] != '/') {
            output += "/";
        }
        // Check output directory exists or not
        if (stat(output.c_str(), &st) == -1) {
#if defined(_WIN32) || defined(_WIN64)
            _mkdir(output.c_str());
#else
            mkdir(output.c_str(), 0755);
#endif
        }
        // Get all files in input directory
        std::cout << "Compressing files in " << input;
        std::cout << " using " << num_threads << " threads" << std::endl;
        std::cout << "Output directory: " << output << std::endl;
        // Parallelize
        namespace gcs = ::google::cloud::storage;
        auto options = google::cloud::Options{}
            .set<gcs::ConnectionPoolSizeOption>(num_threads)
            .set<google::cloud::storage_experimental::HttpVersionOption>("2.0");
        auto client = gcs::Client(options);
        std::vector<std::string> parts = stringSplit(input, "/");
        if (parts.size() == 1) {
            std::cerr << "Invalid gcs URI" << std::endl;
            return 1;
        }
        std::string bucket_name = parts[1];

        // Filter for splitting input into 10 different processes
        //char filter = parts[2][0];
        int num_tar = num_threads / 2;
        mtar_t tarArray[num_tar];
        std::vector<std::string> tarFiles;
        for (int i = 0; i < num_tar; i++) {
            std::string tarFile = output + "AF2_Uniprot_foldcomp." + std::to_string(i) + ".tar";
            tarFiles.push_back(tarFile);
            mtar_open(&tarArray[i], tarFile.c_str(), "w");
        }

        omp_set_num_threads(num_threads);
#pragma omp parallel
        {
#pragma omp single
            // Get object list from gcs bucket
            for (auto&& object_metadata : client.ListObjects(bucket_name, gcs::Projection::NoAcl(), gcs::MaxResults(100000))) {
                std::string obj_name = object_metadata->name();
                // Set zero padding for ID with 4 digits
#pragma omp task firstprivate(obj_name)
                {
                    // Filter for splitting input into 10 different processes
                    // bool skipFilter = filter != '\0' && obj_name.length() >= 9 && obj_name[8] == filter;
                    bool skipFilter = true;
                    bool allowedSuffix = stringEndsWith(".cif", obj_name) || stringEndsWith(".pdb", obj_name);
                    if (skipFilter && allowedSuffix) {
                        auto reader = client.ReadObject(bucket_name, obj_name);
                        if (!reader.status().ok()) {
                            std::cerr << "Could not read object " << obj_name << std::endl;
                        } else {

                            std::string contents{ std::istreambuf_iterator<char>{reader}, {} };
                            CompressedResidue compRes = CompressedResidue();
                            std::string outputFile = output + getFileWithoutExt(obj_name) + ".fcz";
                            compressFromBufferWithoutWriting(compRes, contents, obj_name);
                            //compRes.writeTar(tar, outputFile, compRes.getSize());
                            int tar_id = omp_get_thread_num() / 2;
                            compRes.writeTar(tarArray[tar_id], outputFile, compRes.getSize());
                        }
                    }

                }
            }
            // Close tar
#pragma omp taskwait
            for (int i = 0; i < num_threads; i++) {
                mtar_finalize(&tarArray[i]);
                mtar_close(&tarArray[i]);
            }

        }
#endif
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
