/**
 * File: main.cpp
 * Project: foldcomp
 * Created: 2021-12-23 17:44:53
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Contributor: Milot Mirdita (milot@mirdita.de), Martin Steinegger (themartinsteinegger@gmail.com)
 * Description:
 *     This code contains main function for "foldcomp".
 *     Foldcomp is a fast lossy compression algorithm for protein structure.
 *     It encodes torsion angles with optimal number of bits and reconstruct
 *     3D coordinates from the encoded angles.
 * Usage:
 *    foldcomp compress input.pdb output.fcz
 *    foldcomp decompress input.fcz output.pdb
 * ---
 * Last Modified: 2022-12-05 22:50:47
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
// Headers in the project
#include "atom_coordinate.h"
#include "foldcomp.h"
#include "structure_reader.h"
#include "utility.h"
#include "database_writer.h"
#include "tcbspan.h"
#include "execution_timer.h"

// Standard libraries
#include <cstring>
#include <fstream> // IWYU pragma: keep
#ifdef _WIN32
#include <direct.h>
#include "windows/getopt.h"
#include "windows/dirent.h"
#else
#include <getopt.h>
#endif
#include <iostream>
#include <sstream> // IWYU pragma: keep
#include <string>
#include <vector>

#include <sys/stat.h>

// OpenMP for parallelization
#ifdef OPENMP
#include <omp.h>
#endif

#ifdef HAVE_GCS
#include "google/cloud/storage/client.h"
#endif

static int use_alt_order = 0;
static int anchor_residue_threshold = DEFAULT_ANCHOR_THRESHOLD;
static int save_as_tar = 0;
static int ext_mode = 0;
static int ext_merge = 1;

int print_usage(void) {
    std::cout << "Usage: foldcomp compress <pdb_file> [<fcz_file>]" << std::endl;
    std::cout << "       foldcomp compress [-t number] <pdb_dir|tar> [<fcz_dir>]" << std::endl;
    std::cout << "       foldcomp decompress <fcz_file|tar> [<pdb_file>]" << std::endl;
    std::cout << "       foldcomp decompress [-t number] <fcz_dir|tar> [<pdb_dir>]" << std::endl;
    std::cout << "       foldcomp extract [--plddt|--amino-acid] <fcz_file> [<fasta_file>]" << std::endl;
    std::cout << "       foldcomp extract [--plddt|--amino-acid] [-t number] <fcz_dir|tar> [<fasta_dir>]" << std::endl;
    std::cout << "       foldcomp check <fcz_file>" << std::endl;
    std::cout << "       foldcomp check [-t number] <fcz_dir|tar>" << std::endl;
    std::cout << "       foldcomp rmsd <pdb1|cif1> <pdb2|cif2>" << std::endl;
    std::cout << " -h, --help           print this help message" << std::endl;
    std::cout << " -t, --threads        threads for (de)compression of folders/tar files [default=1]" << std::endl;
    std::cout << " -r, --recursive      recursively look for files in directory [default=0]" << std::endl;
    std::cout << " -f, --file           input is a list of files [default=0]" << std::endl;
    std::cout << " -a, --alt            use alternative atom order [default=false]" << std::endl;
    std::cout << " -b, --break          interval size to save absolute atom coordinates [default=" << anchor_residue_threshold << "]" << std::endl;
    std::cout << " -z, --tar            save as tar file [default=false]" << std::endl;
    std::cout << " -d, --db             save as database [default=false]" << std::endl;
    std::cout << " --skip-discontinuous skip PDB with with discontinuous residues (only batch compression)" << std::endl;
    std::cout << " --plddt              extract pLDDT score (only for extraction mode)" << std::endl;
    std::cout << " --fasta              extract amino acid sequence (only for extraction mode)" << std::endl;
    std::cout << " --no-merge           do not merge output files (only for extraction mode)" << std::endl;
    std::cout << " --time               measure time for compression/decompression" << std::endl;
    return 0;
}

int compress(std::string input, std::string output) {
    StructureReader reader;
    reader.load(input);
    std::vector<AtomCoordinate> atomCoordinates;
    reader.readAllAtoms(atomCoordinates);
    if (atomCoordinates.size() == 0) {
        std::cout << "[Error] No atoms found in the input file: " << input << std::endl;
        return 1;
    }
    std::string title = reader.title;

    // Prototyping for multiple chain support - 2022-10-13 22:01:53
    removeAlternativePosition(atomCoordinates);
    // Identify multiple chains or regions with discontinuous residue indices
    std::vector<std::pair<size_t, size_t>> chain_indices = identifyChains(atomCoordinates);
    // Check if there are multiple chains or regions with discontinuous residue indices
    std::pair<std::string, std::string> outputParts = getFileParts(output);
    std::vector<BackboneChain> compData;
    for (size_t i = 0; i < chain_indices.size(); i++) {
        std::vector<std::pair<size_t, size_t>> frag_indices = identifyDiscontinousResInd(atomCoordinates, chain_indices[i].first, chain_indices[i].second);
        for (size_t j = 0; j < frag_indices.size(); j++) {
            tcb::span<AtomCoordinate> frag_span = tcb::span<AtomCoordinate>(&atomCoordinates[frag_indices[j].first], &atomCoordinates[frag_indices[j].second]);
            Foldcomp compRes;
            compRes.strTitle = title;
            compRes.anchorThreshold = anchor_residue_threshold;
            compData = compRes.compress(frag_span);

            std::string filename;
            if (chain_indices.size() > 1) {
                std::string chain = atomCoordinates[chain_indices[i].first].chain;
                filename = outputParts.first + chain;
            } else {
                filename = outputParts.first;
            }

            if (frag_indices.size() > 1) {
                filename += "_" + std::to_string(j);
            }

            if (outputParts.second != "") {
                filename += "." + outputParts.second;
            }

            // Write compressed data to file
            if (compRes.write(filename) != 0) {
                std::cout << "[Error] Writing file: " << filename << std::endl;
                return -1;
            }
        }
    }
    return 0;
}

int compressWithoutWriting(Foldcomp& compRes, std::string input) {
    StructureReader reader;
    reader.load(input);
    std::vector<AtomCoordinate> atomCoordinates;
    reader.readAllAtoms(atomCoordinates);
    if (atomCoordinates.size() == 0) {
        std::cout << "[Error] No atoms found in the input file: " << input << std::endl;
        return 1;
    }
    std::string title = reader.title;

    removeAlternativePosition(atomCoordinates);

    std::vector<BackboneChain> compData;
    // Convert title to char
    compRes.strTitle = title;
    compRes.anchorThreshold = anchor_residue_threshold;
    compData = compRes.compress(atomCoordinates);
    return 0;
}

int compressFromBufferWithoutWriting(Foldcomp& compRes, const char* data, size_t length, std::string& name) {
    StructureReader reader;
    reader.loadFromBuffer(data, length, name);
    std::vector<AtomCoordinate> atomCoordinates;
    reader.readAllAtoms(atomCoordinates);
    if (atomCoordinates.size() == 0) {
        std::cout << "[Error] No atoms found in the input" << std::endl;
        return 1;
    }
    std::string title = name;

    removeAlternativePosition(atomCoordinates);

    std::vector<BackboneChain> compData;
    // Convert title to char
    compRes.strTitle = name;
    compRes.anchorThreshold = anchor_residue_threshold;
    compData = compRes.compress(atomCoordinates);
    return 0;
}

int decompress(std::istream &file, std::string output) {
    int flag = 0;
    Foldcomp compRes;
    flag = compRes.read(file);
    if (flag != 0) {
        if (flag == -1) {
            std::cerr << "[Error] File is not a valid fcz file" << std::endl;
        } else if (flag == -2) {
            std::cerr << "[Error] Could not restore prevAtoms" << std::endl;
        } else {
            std::cerr << "[Error] Unknown read error" << std::endl;
        }
        return 1;
    }
    std::vector<AtomCoordinate> atomCoordinates;
    compRes.useAltAtomOrder = use_alt_order;
    flag = compRes.decompress(atomCoordinates);
    if (flag != 0) {
        std::cerr << "[Error] decompressing compressed data." << std::endl;
        return 1;
    }
    // Write decompressed data to file
    flag = writeAtomCoordinatesToPDBFile(atomCoordinates, compRes.strTitle, output);
    if (flag != 0) {
        std::cerr << "[Error] Writing decompressed data to file: " << output << std::endl;
        return 1;
    }

    return flag;
}

int extract(std::istream& file, std::string output) {
    int flag = 0;
    Foldcomp compRes;
    flag = compRes.read(file);
    if (flag != 0) {
        if (flag == -1) {
            std::cerr << "[Error] File is not a valid fcz file" << std::endl;
        } else if (flag == -2) {
            std::cerr << "[Error] Could not restore prevAtoms" << std::endl;
        } else {
            std::cerr << "[Error] Unknown read error" << std::endl;
        }
        return 1;
    }
    std::vector<std::string> data;
    compRes.extract(data, ext_mode);
    compRes.writeFASTALike(output, data);
    return 0;
}

int check(std::istream& file, std::string& filename) {
    int flag = 0;
    Foldcomp compRes;
    flag = compRes.read(file);
    if (flag != 0) {
        if (flag == -1) {
            std::cerr << "[Error] File is not a valid fcz file" << std::endl;
        } else if (flag == -2) {
            std::cerr << "[Error] Could not restore prevAtoms" << std::endl;
        } else {
            std::cerr << "[Error] Unknown read error" << std::endl;
        }
        return 1;
    }
    ValidityError err;
    err = compRes.checkValidity();
    printValidityError(err, filename);
    return flag;
}

int rmsd(std::string pdb1, std::string pdb2) {
    // Read
    StructureReader reader1;
    reader1.load(pdb1);
    std::vector<AtomCoordinate> atomCoordinates1;
    reader1.readAllAtoms(atomCoordinates1);
    StructureReader reader2;
    reader2.load(pdb2);
    std::vector<AtomCoordinate> atomCoordinates2;
    reader2.readAllAtoms(atomCoordinates2);
    // Check
    if (atomCoordinates1.size() == 0) {
        std::cerr << "[Error] No atoms found in the input file: " << pdb1 << std::endl;
        return 1;
    }
    if (atomCoordinates2.size() == 0) {
        std::cerr << "[Error] No atoms found in the input file: " << pdb2 << std::endl;
        return 1;
    }
    if (atomCoordinates1.size() != atomCoordinates2.size()) {
        std::cerr << "[Error] The number of atoms in the two files are different." << std::endl;
        return 1;
    }
    std::vector<AtomCoordinate> backbone1 = filterBackbone(atomCoordinates1);
    std::vector<AtomCoordinate> backbone2 = filterBackbone(atomCoordinates2);
    // Print
    std::cout << pdb1 << '\t' << pdb2 << '\t';
    std::cout << backbone1.size() / 3 << '\t' << atomCoordinates1.size() << '\t';
    std::cout << RMSD(backbone1, backbone2) << '\t';
    std::cout << RMSD(atomCoordinates1, atomCoordinates2) << std::endl;
    return 0;
}

int main(int argc, char* const *argv) {
    if (argc < 3) {
        return print_usage();
    }

    int flag = 0;
    int option_index = 0;
    int num_threads = 1;
    int has_output = 0;
    int recursive = 0;
    int file_input = 0;
    int db_output = 0;
    int measure_time = 0;
    int skip_discontinuous = 0;

    // TODO: NEED COMPRESS_MULTIPLE_TAR
    // Mode - non-optional argument
    enum {
        COMPRESS,
        DECOMPRESS,
        COMPRESS_MULTIPLE,
        COMPRESS_MULTIPLE_TAR,
        COMPRESS_MULTIPLE_GCS,
        DECOMPRESS_MULTIPLE,
        DECOMPRESS_MULTIPLE_TAR,
        EXTRACT,
        EXTRACT_MULTIPLE,
        EXTRACT_MULTIPLE_TAR,
        CHECK,
        CHECK_MULTIPLE,
        CHECK_MULTIPLE_TAR,
        RMSD
    } mode = COMPRESS;

    // Define command line options
    static struct option long_options[] = {
            {"help",          no_argument,          0, 'h'},
            {"alt",           no_argument,          0, 'a'},
            {"tar",           no_argument,          0, 'z'},
            {"recursive",     no_argument,          0, 'r'},
            {"file",          no_argument,          0, 'f'},
            {"plddt",         no_argument,  &ext_mode,  0 },
            {"fasta",         no_argument,  &ext_mode,  1 },
            {"no-merge",      no_argument, &ext_merge,  0 },
            {"time",          no_argument, &measure_time, 1 },
            {"skip-discontinuous", no_argument, &skip_discontinuous, 1 },
            {"db",            no_argument,          0, 'd' },
            {"threads", required_argument,          0, 't'},
            {"break",   required_argument,          0, 'b'},
            {0,                         0,          0,  0 }
    };

    // Parse command line options with getopt_long
    flag = getopt_long(argc, argv, "hazrft:b:", long_options, &option_index);

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
            case 'r':
                recursive = 1;
                break;
            case 'f':
                file_input = 1;
                break;
            case 'b':
                anchor_residue_threshold = atoi(optarg);
                break;
            case 'd':
                db_output = 1;
                break;
            case '?':
                return print_usage();
            default:
                break;
        }
        flag = getopt_long(argc, argv, "hazrt:b:", long_options, &option_index);
    }

    // Parse non-option arguments
    // argv[optind]: MODE
    // argv[optind + 1]: INPUT
    // argv[optind + 2]: OUTPUT (optional)

    if ((optind + 1) >= argc) {
        std::cerr << "[Error] Not enough arguments." << std::endl;
        return print_usage();
    }

    struct stat st;
    int fileExists = stat(argv[optind + 1], &st);
    // get mode from command line
    if (strcmp(argv[optind], "compress") == 0) {
        // Check argv[2] is file, directory, or gcs URI
        // TODO: stdin support
        // TODO: COMPRESS_MULTIPLE_TAR
        // If gcs URI, mode = COMPRESS_MULTIPLE_GCS
        // If directory, mode = COMPRESS_MULTIPLE
        // If file, mode = COMPRESS
        char* end = strrchr(argv[optind + 1], '.');
#ifdef HAVE_GCS
        if ((optind + 1) < argc && stringStartsWith("gcs://", argv[optind + 1])) {
            mode = COMPRESS_MULTIPLE_GCS;
            fileExists = 0;
        } else
#endif
        if (file_input || st.st_mode & S_ISDIR(st.st_mode)) {
            mode = COMPRESS_MULTIPLE;
        } else if (end != NULL && strcmp(end, ".tar") == 0) {
            mode = COMPRESS_MULTIPLE_TAR;
        } else {
            mode = COMPRESS;
        }
    } else if (strcmp(argv[optind], "decompress") == 0) {
        // Check argv[2] is file or directory
        // If directory, mode = DECOMPRESS_MULTIPLE
        // If file, mode = DECOMPRESS
        char *end = strrchr(argv[optind + 1], '.');
        if (st.st_mode & S_ISDIR(st.st_mode)) {
            mode = DECOMPRESS_MULTIPLE;
        } else if (end != NULL && strcmp(end, ".tar") == 0) {
            mode = DECOMPRESS_MULTIPLE_TAR;
        } else {
            mode = DECOMPRESS;
        }
    } else if (strcmp(argv[optind], "extract") == 0) {
        // Check argv[2] is file or directory
        // If directory, mode = EXTRACT_MULTIPLE
        // If file, mode = EXTRACT
        char* end = strrchr(argv[optind + 1], '.');
        if (st.st_mode & S_ISDIR(st.st_mode)) {
            mode = EXTRACT_MULTIPLE;
        } else if (end != NULL && strcmp(end, ".tar") == 0) {
            mode = EXTRACT_MULTIPLE_TAR;
        } else {
            mode = EXTRACT;
        }
    } else if (strcmp(argv[optind], "check") == 0){
        char* end = strrchr(argv[optind + 1], '.');
        if (st.st_mode & S_ISDIR(st.st_mode)) {
            mode = CHECK_MULTIPLE;
        }
        else if (end != NULL && strcmp(end, ".tar") == 0) {
            mode = CHECK_MULTIPLE_TAR;
        }
        else {
            mode = CHECK;
        }
    } else if (strcmp(argv[optind], "rmsd") == 0) {
        // NO MULTIPLE MODE FOR RMSD
        mode = RMSD;
    } else {
        return print_usage();
    }
    // Error if no input file given
    if (mode != COMPRESS_MULTIPLE_GCS && fileExists == -1) {
        std::cerr << "[Error] " << argv[optind + 1] << " does not exist." << std::endl;
        return 1;
    }

    std::string input = argv[optind + 1];
    std::string output;
    if (argc == optind + 3) {
        has_output = 1;
        output = argv[optind + 2];
    }

    if (input.back() == '/') {
        input.pop_back();
    }

    if (has_output && output.back() == '/') {
        output.pop_back();
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
        std::ifstream inputFile(input, std::ios::binary);
        // Check if file is open
        if (!inputFile) {
            std::cout << "[Error] Could not open file " << input << std::endl;
            return -1;
        }

        std::cout << "Decompressing " << input << " to " << output << std::endl;
        decompress(inputFile, output);
        inputFile.close();
        flag = 0;
    } else if (mode == EXTRACT) {
        // In extract mode, specific information is directly extracted from the fcz file
        if (!has_output) {
            if (ext_mode == 0) {
                output = getFileWithoutExt(input) + ".plddt.txt";
                std::cout << "Extracting PLDDT from " << input << " to " << output << std::endl;
            } else if (ext_mode == 1) {
                output = getFileWithoutExt(input) + ".fasta";
                std::cout << "Extracting amino acid sequence from " << input << " to " << output << std::endl;
            }
        }
        std::ifstream inputFile(input, std::ios::binary);
        // Check if file is open
        if (!inputFile) {
            std::cout << "[Error] Could not open file " << input << std::endl;
            return -1;
        }
        extract(inputFile, output);
        inputFile.close();
        flag = 0;
    } else if (mode == CHECK){
        // Check if the file is a valid fcz file
        std::ifstream inputFile(input, std::ios::binary);
        std::cout << "Checking " << input << std::endl;
        if (!inputFile) {
            std::cerr << "[Error] Could not open file " << input << std::endl;
            return -1;
        }
        flag = check(inputFile, input);
        inputFile.close();
    } else if (mode == RMSD) {
        // Calculate RMSD between two PDB files
        rmsd(input, output);
        flag = 0;
    } else if (mode == COMPRESS_MULTIPLE || mode == COMPRESS_MULTIPLE_TAR) {
        // compress multiple files
        if (!has_output) {
            if (db_output) {
                output = input + "_db";
            } else if (save_as_tar) {
                output = input + ".fcz.tar";
            } else {
                output = input + "_fcz/";
            }
        } else {
            if (stringEndsWith(".tar", output)) {
                save_as_tar = 1;
            }
        }

        // output variants
        void* handle;
        mtar_t tar_out;
        if (save_as_tar) {
            mtar_open(&tar_out, output.c_str(), "w");
        } else if (db_output) {
            handle = make_writer(output.c_str(), (output + ".index").c_str());
        } else {
            if (stat(output.c_str(), &st) == -1) {
#ifdef _WIN32
                _mkdir(output.c_str());
#else
                mkdir(output.c_str(), 0755);
#endif
            }
        }

        std::cout << "Compressing files in " << input;
        std::cout << " using " << num_threads << " threads" << std::endl;
        if (db_output) {
            std::cout << "Output database: " << output << std::endl;
        } else if (save_as_tar) {
            std::cout << "Output tar file: " << output << std::endl;
        } else {
            std::cout << "Output directory: " << output << std::endl;
        }

        std::vector<std::string> inputs;
        if (file_input) {
            std::ifstream inputFile(input);
            if (!inputFile) {
                std::cerr << "[Error] Could not open file " << input << std::endl;
                return -1;
            }
            std::string line;
            while (std::getline(inputFile, line)) {
                inputs.push_back(line);
            }
        } else {
            inputs.push_back(input);
        }

        for (const std::string& input : inputs) {
            const bool is_tar = stringEndsWith(".tar", input);
            // input variants
            std::vector<std::string> files;
            mtar_t tar_in;
            if (is_tar) {
                if (mtar_open(&tar_in, input.c_str(), "r") != MTAR_ESUCCESS) {
                    if (file_input) {
                        std::cerr << "[Warning] open tar " << input << " failed." << std::endl;
                        continue;
                    } else {
                        std::cerr << "[Error] open tar " << input << " failed." << std::endl;
                        return 1;
                    }
                }
            } else {
                files = getFilesInDirectory(input, recursive);
            }

            unsigned int key = 0;
            if (is_tar) {
#pragma omp parallel shared(tar_in) num_threads(num_threads)
                {
                    std::vector<AtomCoordinate> atomCoordinates;
                    std::vector<BackboneChain> compData;
                    StructureReader reader;

                    bool proceed = true;
                    mtar_header_t header;
                    size_t bufferSize = 1024 * 1024;
                    char* dataBuffer = (char*)malloc(bufferSize);
                    std::string name;
                    while (proceed) {
                        bool writeEntry = true;
#pragma omp critical
                        {
                            if (mtar_read_header(&tar_in, &header) != MTAR_ENULLRECORD) {
                                //TODO GNU tar has special blocks for long filenames
                                name = header.name;
                                if (header.size > bufferSize) {
                                    bufferSize = header.size * 1.5;
                                    dataBuffer = (char*)realloc(dataBuffer, bufferSize);
                                }
                                if (mtar_read_data(&tar_in, dataBuffer, header.size) != MTAR_ESUCCESS) {
                                    std::cerr << "[Error] reading tar entry " << name << " failed." << std::endl;
                                    writeEntry = false;
                                    proceed = false;
                                }
                                else {
                                    writeEntry = true;
                                    proceed = true;
                                }
                                mtar_next(&tar_in);
                                writeEntry = (header.type == MTAR_TREG) ? writeEntry : false;
                            }
                            else {
                                proceed = false;
                                writeEntry = false;
                            }
                        } // end read in
                        if (proceed && writeEntry) {
                            TimerGuard guard(name, measure_time);
                            std::string base = baseName(name);
                            std::pair<std::string, std::string> outputParts = getFileParts(base);
                            std::string outputFile;
                            if (save_as_tar || db_output) {
                                outputFile = outputParts.first;
                            } else if (stringEndsWith("/", output)) {
                                outputFile = output + outputParts.first;
                            } else {
                                outputFile = output + "/" + outputParts.first;
                            }
                            reader.loadFromBuffer(dataBuffer, header.size, outputParts.first);
                            reader.readAllAtoms(atomCoordinates);
                            if (atomCoordinates.size() == 0) {
                                std::cout << "[Error] No atoms found in the input file: " << base << std::endl;
                                continue;
                            }
                            std::string title = reader.title;

                            removeAlternativePosition(atomCoordinates);

                            std::vector<std::pair<size_t, size_t>> chain_indices = identifyChains(atomCoordinates);
                            // Check if there are multiple chains or regions with discontinous residue indices
                            for (size_t i = 0; i < chain_indices.size(); i++) {
                                std::vector<std::pair<size_t, size_t>> frag_indices = identifyDiscontinousResInd(atomCoordinates, chain_indices[i].first, chain_indices[i].second);
                                if (skip_discontinuous && frag_indices.size() > 1) {
                                    std::cout << "[Warning] Skipping discontinuous chain: " << base << std::endl;
                                    continue;
                                }
                                for (size_t j = 0; j < frag_indices.size(); j++) {
                                    tcb::span<AtomCoordinate> frag_span = tcb::span<AtomCoordinate>(&atomCoordinates[frag_indices[j].first], &atomCoordinates[frag_indices[j].second]);
                                    Foldcomp compRes;
                                    compRes.strTitle = title;
                                    compRes.anchorThreshold = anchor_residue_threshold;
                                    compData = compRes.compress(frag_span);

                                    std::string filename;
                                    if (chain_indices.size() > 1) {
                                        std::string chain = atomCoordinates[chain_indices[i].first].chain;
                                        filename = outputFile + chain;
                                    } else {
                                        filename = outputFile;
                                    }

                                    if (frag_indices.size() > 1) {
                                        filename += "_" + std::to_string(j);
                                    }

                                    if (!save_as_tar && !db_output) {
                                        if (outputParts.second != "") {
                                            filename += "." + outputParts.second;
                                        } else {
                                            filename += ".fcz";
                                        }
                                    }

                                    if (db_output) {
                                        std::ostringstream oss;
                                        compRes.writeStream(oss);
#pragma omp critical
                                        {
                                            writer_append(handle, oss.str().c_str(), oss.str().size(), key, filename.c_str());
                                            key++;
                                        }
                                    } else if (save_as_tar) {
#pragma omp critical
                                        {
                                            compRes.writeTar(tar_out, baseName(filename), compRes.getSize());
                                        }
                                    } else {
                                        compRes.write(filename);
                                    }
                                    compData.clear();
                                }
                            }
                            atomCoordinates.clear();
                        }
                    } // end while loop
                    free(dataBuffer);
                }
            } else {
                std::string orig_output = output;
#pragma omp parallel num_threads(num_threads)
                {
                    std::vector<AtomCoordinate> atomCoordinates;
                    std::vector<BackboneChain> compData;
                    StructureReader reader;
#pragma omp for
                    for (size_t i = 0; i < files.size(); i++) {
                        TimerGuard guard(files[i], measure_time);
                        reader.load(files[i]);
                        reader.readAllAtoms(atomCoordinates);
                        if (atomCoordinates.size() == 0) {
                            std::cout << "[Error] No atoms found in the input file: " << files[i] << std::endl;
                            continue;
                        }
                        std::string title = reader.title;

                        removeAlternativePosition(atomCoordinates);

                        std::vector<std::pair<size_t, size_t>> chain_indices = identifyChains(atomCoordinates);
                        // Check if there are multiple chains or regions with discontinous residue indices
                        std::pair<std::string, std::string> outputParts = getFileParts(baseName(files[i]));
                        std::string outputFile;
                        if (save_as_tar || db_output) {
                            outputFile = outputParts.first;
                        } else if (stringEndsWith("/", output)) {
                            outputFile = output + outputParts.first;
                        } else {
                            outputFile = output + "/" + outputParts.first;
                        }
                        for (size_t i = 0; i < chain_indices.size(); i++) {
                            std::vector<std::pair<size_t, size_t>> frag_indices = identifyDiscontinousResInd(atomCoordinates, chain_indices[i].first, chain_indices[i].second);
                            if (skip_discontinuous && frag_indices.size() > 1) {
                                std::cout << "[Warning] Skipping discontinuous chain: " << files[i] << std::endl;
                                continue;
                            }
                            for (size_t j = 0; j < frag_indices.size(); j++) {
                                tcb::span<AtomCoordinate> frag_span = tcb::span<AtomCoordinate>(&atomCoordinates[frag_indices[j].first], &atomCoordinates[frag_indices[j].second]);
                                Foldcomp compRes;
                                compRes.strTitle = title;
                                compRes.anchorThreshold = anchor_residue_threshold;
                                compData = compRes.compress(frag_span);

                                std::string filename;
                                if (chain_indices.size() > 1) {
                                    std::string chain = atomCoordinates[chain_indices[i].first].chain;
                                    filename = outputFile + chain;
                                } else {
                                    filename = outputFile;
                                }

                                if (frag_indices.size() > 1) {
                                    filename += "_" + std::to_string(j);
                                }

                                if (!save_as_tar && !db_output) {
                                    if (outputParts.second != "") {
                                        filename += "." + outputParts.second;
                                    } else {
                                        filename += ".fcz";
                                    }
                                }

                                if (db_output) {
                                    std::ostringstream oss;
                                    compRes.writeStream(oss);
#pragma omp critical
                                    {
                                        writer_append(handle, oss.str().c_str(), oss.str().size(), key, filename.c_str());
                                        key++;
                                    }
                                } else if (save_as_tar) {
#pragma omp critical
                                    {
                                        compRes.writeTar(tar_out, baseName(filename), compRes.getSize());
                                    }
                                } else {
                                    compRes.write(filename);
                                }
                                compData.clear();
                            }
                        }
                        atomCoordinates.clear();
                    }
                }
            }
            if (is_tar) {
                mtar_close(&tar_in);
            }
        }
        if (db_output) {
            free_writer(handle);
        } else if (save_as_tar) {
            mtar_finalize(&tar_out);
            mtar_close(&tar_out);
        }

        flag = 0;
    } else if (mode == COMPRESS_MULTIPLE_GCS) {
        // compress multiple files from gcs
#ifdef HAVE_GCS
        if (!has_output) {
            std::cerr << "Please specify output directory" << std::endl;
            return 1;
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
        int num_tar = num_threads;
        mtar_t tarArray[num_tar];
        std::vector<std::string> tarFiles;
        for (int i = 0; i < num_tar; i++) {
            std::string tarFile = output + "/" + "AF2_Uniprot_foldcomp." + std::to_string(i) + ".tar";
            tarFiles.push_back(tarFile);
            mtar_open(&tarArray[i], tarFile.c_str(), "w");
        }

#pragma omp parallel num_threads(num_threads)
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
                            Foldcomp compRes;
                            std::string outputFile = output + "/" + getFileWithoutExt(obj_name) + ".fcz";
                            compressFromBufferWithoutWriting(compRes, contents.c_str(), contents.length(), obj_name);
                            //compRes.writeTar(tar, outputFile, compRes.getSize());
                            int tar_id = 0;
#ifdef OPENMP
                            tar_id = omp_get_thread_num();
#endif

                            compRes.writeTar(tarArray[tar_id], outputFile, compRes.getSize());
                        }
                    }

                }
            }
            // Close tar
#pragma omp taskwait
            for (int i = 0; i < num_tar; i++) {
                mtar_finalize(&tarArray[i]);
                mtar_close(&tarArray[i]);
            }

        }
#endif
        flag = 0;
    } else if (mode == DECOMPRESS_MULTIPLE || mode == DECOMPRESS_MULTIPLE_TAR) {
        if (!has_output) {
            output = input + "_pdb/";
        }
        // Check output directory exists or not
        if (stat(output.c_str(), &st) == -1) {
#if defined(_WIN32) || defined(_WIN64)
            _mkdir(output.c_str());
#else
            mkdir(output.c_str(), 0755);
#endif
        }
        if (mode == DECOMPRESS_MULTIPLE) {
            // decompress multiple files
            std::cout << "Decompressing files in " << input;
            std::cout << " using " << num_threads << " threads" << std::endl;
            std::cout << "Output directory: " << output << std::endl;
            std::vector<std::string> files = getFilesInDirectory(input, recursive);
#pragma omp parallel num_threads(num_threads)
            {
#pragma omp for
                for (size_t i = 0; i < files.size(); i++) {
                    TimerGuard guard(files[i], measure_time);
                    std::ifstream input(files[i], std::ios::binary);
                    // Check if file is open
                    if (!input) {
                        std::cout << "[Error] Could not open file " << files[i] << std::endl;
                        continue;
                    }
                    std::string outputFile = output + "/" + baseName(getFileWithoutExt(files[i])) + ".pdb";
                    decompress(input, outputFile);
                }
            }
            flag = 0;
        } else if (mode == DECOMPRESS_MULTIPLE_TAR) {
            mtar_t tar;
            if (mtar_open(&tar, input.c_str(), "r") != MTAR_ESUCCESS) {
                std::cerr << "[Error] open tar " << input.c_str() << " failed." << std::endl;
                return 1;
            }
#pragma omp parallel shared(tar) num_threads(num_threads)
            {
                bool proceed = true;
                mtar_header_t header;
                size_t bufferSize = 1024 * 1024;
                char *dataBuffer = (char *) malloc(bufferSize);
                std::string name;
                while (proceed) {
                    bool writeEntry = true;
#pragma omp critical
                    {
                        if (mtar_read_header(&tar, &header) != MTAR_ENULLRECORD) {
                            //TODO GNU tar has special blocks for long filenames
                            name = header.name;
                            if (header.size > bufferSize) {
                                bufferSize = header.size * 1.5;
                                dataBuffer = (char *) realloc(dataBuffer, bufferSize);
                            }
                            if (mtar_read_data(&tar, dataBuffer, header.size) != MTAR_ESUCCESS) {
                                std::cerr << "[Error] reading tar entry " << name << " failed." << std::endl;
                                writeEntry = false;
                                proceed = false;
                            } else {
                                writeEntry = true;
                                proceed = true;
                            }
                            mtar_next(&tar);
                            writeEntry = (header.type == MTAR_TREG) ? writeEntry : false;
                        } else {
                            proceed = false;
                            writeEntry = false;
                        }
                    } // end read in
                    if (proceed && writeEntry) {
                        std::istringstream input(std::string(dataBuffer, header.size));
                        std::string name_clean = baseName(name);
                        std::string outputFile = output + "/" + name_clean + ".pdb";
                        decompress(input, outputFile);
                    }
                } // end while loop
            } // end openmp
            mtar_close(&tar);
            flag = 0;
        }
    } else if (mode == EXTRACT_MULTIPLE || mode == EXTRACT_MULTIPLE_TAR) {
            // extract multiple files
            if (mode == EXTRACT_MULTIPLE) {
                if (!has_output) {
                    if (ext_mode == 0) {
                        output = input + "_plddt/";
                    } else if (ext_mode == 1) {
                        output = input + "_fasta/";
                    }
                }
                // Check output directory exists or not
                if (stat(output.c_str(), &st) == -1) {
#if defined(_WIN32) || defined(_WIN64)
                    _mkdir(output.c_str());
#else
                    mkdir(output.c_str(), 0755);
#endif
                }
                std::string defaultOutputFile = "";
                if (ext_mode == 0){
                    defaultOutputFile = output + "plddt.txt";
                } else if (ext_mode == 1) {
                    defaultOutputFile = output + "aa.fasta";
                }
                std::ofstream defaultOutput;
                if (ext_merge == 1) {
                    defaultOutput.open(defaultOutputFile, std::ios::out);
                }
                std::cout << "Extracting files in " << input;
                std::cout << " using " << num_threads << " threads" << std::endl;
                std::cout << "Output directory: " << output << std::endl;
                std::vector<std::string> files = getFilesInDirectory(input, recursive);
#pragma omp parallel num_threads(num_threads)
                {
                    std::string buffer;
                    buffer.reserve(1024 * 1024);
                    std::vector<std::string> data;
                    data.reserve(1024);
#pragma omp for
                    for (size_t i = 0; i < files.size(); i++) {
                        std::ifstream input(files[i], std::ios::binary);
                        if (!input) {
                            std::cout << "[Error] Could not open file " << files[i] << std::endl;
                            continue;
                        }
                        if (ext_merge == 1) {
                            Foldcomp compRes;
                            int flag = compRes.read(input);
                            if (flag != 0) {
                                if (flag == -1) {
                                    std::cerr << "[Error] File is not a valid fcz file" << std::endl;
                                } else if (flag == -2) {
                                    std::cerr << "[Error] Could not restore prevAtoms" << std::endl;
                                } else {
                                    std::cerr << "[Error] Unknown read error" << std::endl;
                                }
                                continue;
                            }
                            compRes.extract(data, ext_mode);
                            buffer.append(1, '>');
                            buffer.append(compRes.strTitle);
                            buffer.append(1, '\n');
                            for (size_t j = 0; j < data.size(); j++) {
                                buffer.append(data[j]);
                            }
                            buffer.append(1, '\n');
                        } else {
                            std::string outputFile;
                            if (ext_mode == 0) {
                                // output file extension is ".plddt.txt"
                                outputFile = output + "/" + baseName(getFileWithoutExt(files[i])) + ".plddt.txt";
                            }
                            else if (ext_mode == 1) {
                                outputFile = output + "/" + baseName(getFileWithoutExt(files[i])) + ".fasta";
                            }
                            extract(input, outputFile);
                        }
                        input.close();
                        data.clear();
                        if (ext_merge == 1 && buffer.size() > 1024 * 1024) {
                            #pragma omp critical
                            {
                                defaultOutput << buffer;
                            }
                            buffer.clear();
                        }
                    }
                    if (ext_merge == 1 && buffer.size() > 0) {
                        #pragma omp critical
                        {
                            defaultOutput << buffer;
                        }
                        buffer.clear();
                    }
                }
                flag = 0;
                if (ext_merge == 1) {
                    // close and delete defaultOutput
                    defaultOutput.close();
                }
            } else if (mode == EXTRACT_MULTIPLE_TAR) {
                mtar_t tar;
                if (mtar_open(&tar, input.c_str(), "r") != MTAR_ESUCCESS) {
                    std::cerr << "[Error] Open tar " << input << " failed." << std::endl;
                    return 1;
                }
                std::string defaultOutputFile = "";
                if (!has_output) {
                    if (ext_mode == 0) {
                        output = getFileWithoutExt(input) + "_plddt/";
                        defaultOutputFile = output + "plddt.txt";
                    }
                    else if (ext_mode == 1) {
                        output = getFileWithoutExt(input) + "_fasta/";
                        defaultOutputFile = output + "aa.fasta";
                    }
                }

                // Check output directory exists or not
                if (stat(output.c_str(), &st) == -1) {
#if defined(_WIN32) || defined(_WIN64)
                    _mkdir(output.c_str());
#else
                    mkdir(output.c_str(), 0755);
#endif
                }
                std::ofstream defaultOutput;
                if (ext_merge == 1) {
                    defaultOutput.open(defaultOutputFile, std::ios::out);
                }
                std::cout << "Extracting files in " << input << " using " << num_threads << " threads" << std::endl;
                // TAR READING PART BY MARTIN STEINEGGER
#pragma omp parallel shared(tar) num_threads(num_threads)
                {
                    bool proceed = true;
                    mtar_header_t header;
                    size_t bufferSize = 1024 * 1024;
                    char* dataBuffer = (char*)malloc(bufferSize);
                    std::string buffer;
                    buffer.reserve(1024 * 1024);
                    std::string name;
                    while (proceed) {
                        bool writeEntry = true;
#pragma omp critical
                        {
                            if (mtar_read_header(&tar, &header) != MTAR_ENULLRECORD) {
                                //TODO GNU tar has special blocks for long filenames
                                name = header.name;
                                if (header.size > bufferSize) {
                                    bufferSize = header.size * 1.5;
                                    dataBuffer = (char*)realloc(dataBuffer, bufferSize);
                                }
                                if (mtar_read_data(&tar, dataBuffer, header.size) != MTAR_ESUCCESS) {
                                    std::cerr << "[Error] Reading tar entry " << name << " failed." << std::endl;
                                    writeEntry = false;
                                    proceed = false;
                                }
                                else {
                                    writeEntry = true;
                                    proceed = true;
                                }
                                mtar_next(&tar);
                                writeEntry = (header.type == MTAR_TREG) ? writeEntry : false;
                            }
                            else {
                                proceed = false;
                                writeEntry = false;
                            }
                        } // end read in
                        if (proceed && writeEntry) {
                            std::istringstream input(std::string(dataBuffer, header.size));
                            if (ext_merge == 1) {
                                std::vector<std::string> data;
                                Foldcomp compRes;
                                int flag = compRes.read(input);
                                if (flag != 0) {
                                    if (flag == -1) {
                                        std::cerr << "[Error] File is not a valid fcz file" << std::endl;
                                    } else if (flag == -2) {
                                        std::cerr << "[Error] Could not restore prevAtoms" << std::endl;
                                    } else {
                                        std::cerr << "[Error] Unknown read error" << std::endl;
                                    }
                                    continue;
                                }
                                compRes.extract(data, ext_mode);
                                buffer.append(1, '>');
                                buffer.append(compRes.strTitle);
                                buffer.append(1, '\n');
                                for (size_t j = 0; j < data.size(); j++) {
                                    buffer.append(data[j]);
                                }
                                buffer.append(1, '\n');
                                #pragma omp critical
                                {
                                    defaultOutput << buffer;
                                }
                                buffer.clear();
                            } else {
                                std::string outputFile;
                                std::string name_clean = baseName(name);
                                if (ext_mode == 0) {
                                    // output file extension is ".plddt.txt"
                                    outputFile = output + "/" + name_clean + ".plddt.txt";
                                }
                                else if (ext_mode == 1) {
                                    outputFile = output + "/" + name_clean + ".fasta";
                                }
                                extract(input, outputFile);
                            }
                        }
                    } // end while loop
                    free(dataBuffer);
                } // end openmp
                mtar_close(&tar);
                flag = 0;
                if (ext_merge == 1) {
                    defaultOutput.close();
                }
            }
    } else if (mode == CHECK_MULTIPLE) {
        std::cout << "Checking files in " << input << " using " << num_threads << " threads" << std::endl;
        std::vector<std::string> files = getFilesInDirectory(input, recursive);
#pragma omp parallel num_threads(num_threads)
        {
#pragma omp for
            for (size_t i = 0; i < files.size(); i++) {
                std::ifstream input(files[i], std::ios::binary);
                // Check if file is open
                if (!input) {
                    std::cerr << "[Error] Could not open file " << files[i] << std::endl;
                    continue;
                }
                check(input, files[i]);
                input.close();
            }
        }
        flag = 0;
    } else if (mode == CHECK_MULTIPLE_TAR) {
        mtar_t tar;
        if (mtar_open(&tar, input.c_str(), "r") != MTAR_ESUCCESS) {
            std::cerr << "[Error] open tar " << input << " failed." << std::endl;
            return 1;
        }
        std::cout << "Checking files in " << input << " using " << num_threads << " threads" << std::endl;
        // TAR READING PART BY MARTIN STEINEGGER
#pragma omp parallel shared(tar) num_threads(num_threads)
        {
            bool proceed = true;
            mtar_header_t header;
            size_t bufferSize = 1024 * 1024;
            char* dataBuffer = (char*)malloc(bufferSize);
            std::string name;
            while (proceed) {
                bool writeEntry = true;
#pragma omp critical
                {
                    if (mtar_read_header(&tar, &header) != MTAR_ENULLRECORD) {
                        //TODO GNU tar has special blocks for long filenames
                        name = header.name;
                        if (header.size > bufferSize) {
                            bufferSize = header.size * 1.5;
                            dataBuffer = (char*)realloc(dataBuffer, bufferSize);
                        }
                        if (mtar_read_data(&tar, dataBuffer, header.size) != MTAR_ESUCCESS) {
                            std::cerr << "[Error] Reading tar entry " << name << " failed." << std::endl;
                            writeEntry = false;
                            proceed = false;
                        }
                        else {
                            writeEntry = true;
                            proceed = true;
                        }
                        mtar_next(&tar);
                        writeEntry = (header.type == MTAR_TREG) ? writeEntry : false;
                    }
                    else {
                        proceed = false;
                        writeEntry = false;
                    }
                } // end read in
                if (proceed && writeEntry) {
                    std::istringstream input(std::string(dataBuffer, header.size));
                    std::string name_clean = baseName(name);
                    check(input, name_clean);
                }
            } // end while loop
        } // end openmp
        flag = 0;
    } else {
        std::cerr << "Invalid mode." << std::endl;
        return 1;
    }    // Print log
    if (mode != RMSD) {
        std::cout << "Done." << std::endl;
    }
    return flag;
}
