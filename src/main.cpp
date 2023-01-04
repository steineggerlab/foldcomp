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
#include "input_processor.h"

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

int rmsd(const std::string& pdb1, const std::string& pdb2) {
    StructureReader reader;
    reader.load(pdb1);
    std::vector<AtomCoordinate> atomCoordinates1;
    reader.readAllAtoms(atomCoordinates1);
    reader.load(pdb2);
    std::vector<AtomCoordinate> atomCoordinates2;
    reader.readAllAtoms(atomCoordinates2);
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

    int option_index = 0;
    int num_threads = 1;
    int has_output = 0;
    int recursive = 0;
    int file_input = 0;
    int db_output = 0;
    int measure_time = 0;
    int skip_discontinuous = 0;

    // Mode - non-optional argument
    enum {
        COMPRESS,
        DECOMPRESS,
        EXTRACT,
        CHECK,
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
    int flag = getopt_long(argc, argv, "hazrft:b:", long_options, &option_index);
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

    struct stat inputStat;
    int inputExists = stat(argv[optind + 1], &inputStat);
    const char* outputSuffix = "";
    bool mayHaveOutput = false;
    // get mode from command line
    if (strcmp(argv[optind], "compress") == 0) {
        mode = COMPRESS;
#ifdef HAVE_GCS
        if ((optind + 1) < argc && stringStartsWith("gcs://", argv[optind + 1])) {
            fileExists = 1;
        }
#endif
        mayHaveOutput = true;
        outputSuffix = "fcz";
    } else if (strcmp(argv[optind], "decompress") == 0) {
        mode = DECOMPRESS;
        mayHaveOutput = true;
        outputSuffix = "pdb";
    } else if (strcmp(argv[optind], "extract") == 0) {
        mode = EXTRACT;
        mayHaveOutput = true;
        if (ext_mode == 0){
            outputSuffix = "plddt.txt";
        } else if (ext_mode == 1) {
            outputSuffix = "aa.fasta";
        }
    } else if (strcmp(argv[optind], "check") == 0){
        mode = CHECK;
    } else if (strcmp(argv[optind], "rmsd") == 0) {
        mode = RMSD;
    } else {
        return print_usage();
    }

    // Error if no input file given
    if (inputExists == -1) {
        std::cerr << "[Error] " << argv[optind + 1] << " does not exist." << std::endl;
        return EXIT_FAILURE;
    }

    std::string input = argv[optind + 1];
    while (input.back() == '/') {
        input.pop_back();
    }
    std::vector<std::string> inputs;
    if (file_input) {
        std::ifstream inputFile(input);
        if (!inputFile) {
            std::cerr << "[Error] Could not open file " << input << std::endl;
            return EXIT_FAILURE;
        }
        std::string line;
        while (std::getline(inputFile, line)) {
            inputs.push_back(line);
        }
    } else {
        inputs.push_back(input);
    }

    std::string output;
    if (argc == optind + 3) {
        has_output = 1;
        output = argv[optind + 2];
        if (stringEndsWith(".tar", output)) {
            save_as_tar = 1;
        }
    }

    bool isSingleFileInput = !file_input && (S_ISREG(inputStat.st_mode) || S_ISLNK(inputStat.st_mode));
    if (!file_input) {
        struct stat st;
        if (stringEndsWith(".tar", input) || stringEndsWith(".tar.gz", input) || stringEndsWith(".tgz", input)) {
            isSingleFileInput = false;
        } else if (stat((input + ".dbtype").c_str(), &st) == 0) {
            isSingleFileInput = false;
        }
#ifdef HAVE_GCS
        else if (stringStartsWith("gcs://", input)) {
            isSingleFileInput = false;
        }
#endif
    }

    while (has_output && output.back() == '/') {
        output.pop_back();
    }
    if (mayHaveOutput && !has_output) {
        if (db_output) {
            output = input + "_db";
        } else if (save_as_tar) {
            output = input + "." +  outputSuffix + ".tar";
        } else {
            if (isSingleFileInput) {
                output = getFileWithoutExt(input) + "." + outputSuffix;
            } else {
                // TODO EXTRACT split at .
                output = input + "_" + outputSuffix + "/";
            }
        }
    }

    if (mode == RMSD) {
        // Calculate RMSD between two PDB files
        rmsd(input, output);
    } else if (mode == COMPRESS) {
        // output variants
        void* handle;
        mtar_t tar_out;
        if (save_as_tar) {
            mtar_open(&tar_out, output.c_str(), "w");
        } else if (db_output) {
            handle = make_writer(output.c_str(), (output + ".index").c_str());
        } else if (!isSingleFileInput) {
            struct stat st;
            if (stat(output.c_str(), &st) == -1) {
#ifdef _WIN32
                _mkdir(output.c_str());
#else
                mkdir(output.c_str(), 0755);
#endif
            }
        }

        if (isSingleFileInput) {
            std::cout << "Compressing " << input << " to " << output << std::endl;
        } else {
            std::cout << "Compressing files in " << input;
            std::cout << " using " << num_threads << " threads" << std::endl;
            if (db_output) {
                std::cout << "Output database: " << output << std::endl;
            } else if (save_as_tar) {
                std::cout << "Output tar file: " << output << std::endl;
            } else {
                std::cout << "Output directory: " << output << std::endl;
            }
        }

        for (const std::string& input : inputs) {
            Processor* processor;
            struct stat st;
            if (stringEndsWith(".tar", input) || stringEndsWith(".tar.gz", input) || stringEndsWith(".tgz", input)) {
                processor = new TarProcessor(input);
            } else if (stat((input + ".dbtype").c_str(), &st) == 0) {
                processor = new DatabaseProcessor(input);
            } 
#ifdef HAVE_GCS
            else if (stringStartsWith("gcs://", input)) {
                processor = new GcsProcessor(input);
            }
#endif 
            else {
                processor = new DirectoryProcessor(input, recursive);
            }
            unsigned int key = 0;
            process_entry_func func = [&](const char* name, const char* dataBuffer, size_t size) -> bool {
                TimerGuard guard(name, measure_time);
                std::vector<AtomCoordinate> atomCoordinates;
                std::vector<BackboneChain> compData;
                StructureReader reader;

                std::string base = baseName(name);
                std::pair<std::string, std::string> outputParts = getFileParts(base);
                std::string outputFile;
                if (save_as_tar || db_output) {
                    outputFile = outputParts.first;
                } else if (isSingleFileInput) {
                    outputParts = getFileParts(output);
                    outputFile = outputParts.first;
                } else {
                    outputFile = output + "/" + outputParts.first;
                }
                reader.loadFromBuffer(dataBuffer, size, base);
                reader.readAllAtoms(atomCoordinates);
                if (atomCoordinates.size() == 0) {
                    std::cerr << "[Error] No atoms found in the input file: " << base << std::endl;
                    return false;
                }
                std::string title = outputParts.first;

                removeAlternativePosition(atomCoordinates);

                std::vector<std::pair<size_t, size_t>> chain_indices = identifyChains(atomCoordinates);
                // Check if there are multiple chains or regions with discontinous residue indices
                for (size_t i = 0; i < chain_indices.size(); i++) {
                    std::vector<std::pair<size_t, size_t>> frag_indices = identifyDiscontinousResInd(atomCoordinates, chain_indices[i].first, chain_indices[i].second);
                    if (skip_discontinuous && frag_indices.size() > 1) {
                        std::cerr << "[Warning] Skipping discontinuous chain: " << base << std::endl;
                        continue;
                    }
                    for (size_t j = 0; j < frag_indices.size(); j++) {
                        tcb::span<AtomCoordinate> frag_span = tcb::span<AtomCoordinate>(&atomCoordinates[frag_indices[j].first], atomCoordinates.data() + frag_indices[j].second);
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
                            if (outputParts.second != "" && outputParts.second != "pdb" && outputParts.second != "cif") {
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
                return true;
            };
            processor->run(func, num_threads);
            delete processor;
        }
        if (db_output) {
            free_writer(handle);
        } else if (save_as_tar) {
            mtar_write_finalize(&tar_out);
            mtar_close(&tar_out);
        }
    } else if (mode == DECOMPRESS) {
        void* handle;
        mtar_t tar_out;
        if (save_as_tar) {
            mtar_open(&tar_out, output.c_str(), "w");
        } else if (db_output) {
            handle = make_writer(output.c_str(), (output + ".index").c_str());
        } else if (!isSingleFileInput) {
            struct stat st;
            if (stat(output.c_str(), &st) == -1) {
#ifdef _WIN32
                _mkdir(output.c_str());
#else
                mkdir(output.c_str(), 0755);
#endif
            }
        }

        if (isSingleFileInput) {
            std::cout << "Decompressing " << input << " to " << output << std::endl;
        } else {
            std::cout << "Decompressing files in " << input;
            std::cout << " using " << num_threads << " threads" << std::endl;
            if (db_output) {
                std::cout << "Output database: " << output << std::endl;
            } else if (save_as_tar) {
                std::cout << "Output tar file: " << output << std::endl;
            } else {
                std::cout << "Output directory: " << output << std::endl;
            }
        }

        for (const std::string& input : inputs) {
            Processor* processor;
            struct stat st;
            if (stringEndsWith(".tar", input) || stringEndsWith(".tar.gz", input) || stringEndsWith(".tgz", input)) {
                processor = new TarProcessor(input);
            } else if (stat((input + ".dbtype").c_str(), &st) == 0) {
                processor = new DatabaseProcessor(input);
            } 
#ifdef HAVE_GCS
            else if (stringStartsWith("gcs://", input)) {
                processor = new GcsProcessor(input);
            }
#endif 
            else {
                processor = new DirectoryProcessor(input, recursive);
            }
            unsigned int key = 0;
            process_entry_func func = [&](const char* name, const char* dataBuffer, size_t size) -> bool {
                TimerGuard guard(name, measure_time);
                Foldcomp compRes;
                std::istringstream input(std::string(dataBuffer, size));
                int flag = compRes.read(input);
                if (flag != 0) {
                    if (flag == -1) {
                        std::cerr << "[Error] File is not a valid fcz file" << std::endl;
                    } else if (flag == -2) {
                        std::cerr << "[Error] Could not restore prevAtoms" << std::endl;
                    } else {
                        std::cerr << "[Error] Unknown read error" << std::endl;
                    }
                    return false;
                }
                std::vector<AtomCoordinate> atomCoordinates;
                compRes.useAltAtomOrder = use_alt_order;
                flag = compRes.decompress(atomCoordinates);
                if (flag != 0) {
                    std::cerr << "[Error] decompressing compressed data." << std::endl;
                    return false;
                }

                std::string base = baseName(name);
                std::pair<std::string, std::string> outputParts = getFileParts(base);
                std::string outputFile;
                if (save_as_tar || db_output) {
                    outputFile = outputParts.first;
                } else if (isSingleFileInput) {
                    outputFile = output;
                } else {
                    outputFile = output + "/" + outputParts.first + "." + outputSuffix;
                }

                if (db_output) {
                    std::ostringstream oss;
                    writeAtomCoordinatesToPDB(atomCoordinates, compRes.strTitle, oss);
                    oss << '\0';
#pragma omp critical
                    {
                        writer_append(handle, oss.str().c_str(), oss.str().size(), key, outputFile.c_str());
                        key++;
                    }
                } else if (save_as_tar) {
                    std::ostringstream oss;
                    writeAtomCoordinatesToPDB(atomCoordinates, compRes.strTitle, oss);
#pragma omp critical
                    {
                        std::string os = oss.str();
                        mtar_write_file_header(&tar_out, outputFile.c_str(), os.size());
                        mtar_write_data(&tar_out, os.c_str(), os.size());
                    }
                } else {
                    // Write decompressed data to file
                    flag = writeAtomCoordinatesToPDBFile(atomCoordinates, compRes.strTitle, outputFile);
                    if (flag != 0) {
                        std::cerr << "[Error] Writing decompressed data to file: " << output << std::endl;
                        return false;
                    }
                }
                return true;
            };
            processor->run(func, num_threads);
            delete processor;
        }
        if (db_output) {
            free_writer(handle);
        } else if (save_as_tar) {
            mtar_write_finalize(&tar_out);
            mtar_close(&tar_out);
        }
    } else if (mode == EXTRACT) {
        void* handle;
        mtar_t tar_out;
        if (save_as_tar) {
            mtar_open(&tar_out, output.c_str(), "w");
        } else if (db_output) {
            handle = make_writer(output.c_str(), (output + ".index").c_str());
        } else {
            struct stat st;
            if (stat(output.c_str(), &st) == -1) {
#ifdef _WIN32
                _mkdir(output.c_str());
#else
                mkdir(output.c_str(), 0755);
#endif
            }
        }

        if (isSingleFileInput) {
            std::cout << "Extracting " << input << " to " << output << std::endl;
        } else {
            std::cout << "Extracting files in " << input;
            std::cout << " using " << num_threads << " threads" << std::endl;
            if (db_output) {
                std::cout << "Output database: " << output << std::endl;
            } else if (save_as_tar) {
                std::cout << "Output tar file: " << output << std::endl;
            } else {
                std::cout << "Output directory: " << output << std::endl;
            }
        }

        bool isMergedOutput = false;
        std::ofstream default_out;
        if (!save_as_tar && !db_output && !isSingleFileInput && ext_merge == 0) {
            default_out.open(output, std::ios::out);
            isMergedOutput = true;
        }

        for (const std::string& input : inputs) {
            Processor* processor;
            struct stat st;
            if (stringEndsWith(".tar", input) || stringEndsWith(".tar.gz", input) || stringEndsWith(".tgz", input)) {
                processor = new TarProcessor(input);
            } else if (stat((input + ".dbtype").c_str(), &st) == 0) {
                processor = new DatabaseProcessor(input);
            } 
#ifdef HAVE_GCS
            else if (stringStartsWith("gcs://", input)) {
                processor = new GcsProcessor(input);
            }
#endif 
            else {
                processor = new DirectoryProcessor(input, recursive);
            }
            unsigned int key = 0;
            process_entry_func func = [&](const char* name, const char* dataBuffer, size_t size) -> bool {
                std::istringstream input(std::string(dataBuffer, size));
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
                    return false;
                }
                std::string data;
                compRes.extract(data, ext_mode);
                std::string base = baseName(name);
                std::pair<std::string, std::string> outputParts = getFileParts(base);
                std::string outputFile;
                if (save_as_tar || db_output) {
                    outputFile = outputParts.first;
                } else if (isSingleFileInput) {
                    outputFile = output;
                } else {
                    outputFile = output + "/" + outputParts.first + "." + outputSuffix;
                }

                if (isMergedOutput) {
#pragma omp critical
                    {
                        compRes.writeFASTALike(default_out, data);
                    }
                } else if (db_output) {
                    std::ostringstream oss;
                    compRes.writeFASTALike(oss, data);
                    oss << '\0';
#pragma omp critical
                    {
                        writer_append(handle, oss.str().c_str(), oss.str().size(), key, outputFile.c_str());
                        key++;
                    }
                } else if (save_as_tar) {
                    std::ostringstream oss;
                    compRes.writeFASTALike(oss, data);
#pragma omp critical
                    {
                        std::string os = oss.str();
                        mtar_write_file_header(&tar_out, outputFile.c_str(), os.size());
                        mtar_write_data(&tar_out, os.c_str(), os.size());
                    }
                } else {
                    std::ofstream output(outputFile);
                    if (!output) {
                        std::cerr << "[Error] Could not open file " << outputFile << std::endl;
                        return false;
                    }
                    compRes.writeFASTALike(output, data);
                }

                return true;
            };
            processor->run(func, num_threads);
            delete processor;
        }
        if (db_output) {
            free_writer(handle);
        } else if (save_as_tar) {
            mtar_write_finalize(&tar_out);
            mtar_close(&tar_out);
        }
    } else if (mode == CHECK) {
        if (inputs.size() == 1) {
            std::cout << "Checking " << input << std::endl;
        } else {
            std::cout << "Checking files in " << input;
            std::cout << " using " << num_threads << " threads" << std::endl;
        }

        for (const std::string& input : inputs) {
            Processor* processor;
            struct stat st;
            if (stringEndsWith(".tar", input) || stringEndsWith(".tar.gz", input) || stringEndsWith(".tgz", input)) {
                processor = new TarProcessor(input);
            } else if (stat((input + ".dbtype").c_str(), &st) == 0) {
                processor = new DatabaseProcessor(input);
            } 
#ifdef HAVE_GCS
            else if (stringStartsWith("gcs://", input)) {
                processor = new GcsProcessor(input);
            }
#endif 
            else {
                processor = new DirectoryProcessor(input, recursive);
            }
            process_entry_func func = [&](const char* name, const char* dataBuffer, size_t size) -> bool {
                std::istringstream input(std::string(dataBuffer, size));
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
                    return false;
                }
                ValidityError err = compRes.checkValidity();
                std::string sname(name);
                printValidityError(err, sname);
                return true;
            };
            processor->run(func, num_threads);
            delete processor;
        }
    }
    return EXIT_SUCCESS;
}
