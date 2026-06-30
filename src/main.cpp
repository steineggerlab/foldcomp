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
 * Last Modified: 2026-02-10 20:09:37
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
// Headers in the project
#include "amino_acid.h"
#include "atom_coordinate.h"
#include "foldcomp.h"
#include "structure_reader.h"
#include "utility.h"
#include "database_writer.h"
#include "tcbspan.h"
#include "execution_timer.h"
#include "input_processor.h"
#include "structure_codec.h"

// Standard libraries
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream> // IWYU pragma: keep
#include <iomanip>
#include <limits>
#ifdef _WIN32
#include <direct.h>
#include "windows/getopt.h"
#include "windows/dirent.h"
#else
#include <getopt.h>
#endif
#include <iostream>
#include <sstream> // IWYU pragma: keep
#include <set>
#include <map>
#include <string>
#include <tuple>
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
static int ext_plddt_digits = 1;
static int ext_merge = 1;
static int ext_use_title = 0;
static int overwrite = 0;
static float max_backbone_rmsd = std::numeric_limits<float>::infinity();

static bool isMMCIFOutputPath(const std::string& path) {
    return stringEndsWith(".cif", path) ||
           stringEndsWith(".mmcif", path) ||
           stringEndsWith(".cif.tar", path) ||
           stringEndsWith(".mmcif.tar", path);
}

// version
#define FOLDCOMP_VERSION "1.0.0"

static bool residueNeedsRawFallback(const tcb::span<const AtomCoordinate>& residueAtoms) {
    if (residueAtoms.empty()) {
        return false;
    }
    auto aaIt = Foldcomp::AAS.find(residueAtoms[0].residue);
    if (aaIt == Foldcomp::AAS.end()) {
        return true;
    }

    int matchedCanonicalAtoms = 0;
    bool hasOxt = false;
    for (const auto& atom : residueAtoms) {
        if (atom.altloc != ' ' && atom.altloc != '\0') {
            return true;
        }
        if (atom.insertion_code != ' ' && atom.insertion_code != '\0') {
            return true;
        }
        if (atom.atom == "OXT") {
            if (hasOxt) {
                return true;
            }
            hasOxt = true;
            continue;
        }
        bool found = false;
        for (const auto& canonicalAtom : aaIt->second.atoms) {
            if (atom.atom == canonicalAtom) {
                matchedCanonicalAtoms++;
                found = true;
                break;
            }
        }
        if (!found) {
            return true;
        }
    }
    int expectedAtomCount = static_cast<int>(aaIt->second.atoms.size()) + (hasOxt ? 1 : 0);
    return static_cast<int>(residueAtoms.size()) != expectedAtomCount ||
           matchedCanonicalAtoms != static_cast<int>(aaIt->second.atoms.size());
}

static bool regionNeedsRawFallback(const tcb::span<AtomCoordinate>& atoms) {
    size_t residuesWithOxt = 0;
    size_t residueStart = 0;
    size_t residueCount = 0;
    for (size_t i = 1; i <= atoms.size(); i++) {
        bool endOfResidue = (i == atoms.size()) || startsNewResidue(atoms[i], atoms[i - 1]);
        if (!endOfResidue) {
            continue;
        }
        residueCount++;
        bool hasOxt = false;
        for (size_t j = residueStart; j < i; j++) {
            const auto& atom = atoms[j];
            if (atom.atom == "OXT") {
                hasOxt = true;
                break;
            }
        }
        if (!hasOxt) {
        } else {
            residuesWithOxt++;
            if (residuesWithOxt > 1 || i != atoms.size()) {
                return true;
            }
        }
        if (residueNeedsRawFallback(
                tcb::span<const AtomCoordinate>(atoms.data() + residueStart, i - residueStart))) {
            return true;
        }
        residueStart = i;
    }
    return false;
}

static size_t countResidues(const tcb::span<AtomCoordinate>& atoms) {
    if (atoms.empty()) {
        return 0;
    }
    size_t residueCount = 1;
    for (size_t i = 1; i < atoms.size(); i++) {
        if (startsNewResidue(atoms[i], atoms[i - 1])) {
            residueCount++;
        }
    }
    return residueCount;
}

static bool shouldStoreSmallMixedFragmentAsRaw(
    const tcb::span<AtomCoordinate>& chainSpan,
    const std::vector<BackboneRegion>& regions,
    const std::vector<bool>& regionNeedsRaw
) {
    if (regions.size() <= 1) {
        return false;
    }
    if (countResidues(chainSpan) > 24) {
        return false;
    }
    for (size_t i = 0; i < regions.size(); i++) {
        if (!regions[i].encodable || regionNeedsRaw[i]) {
            return true;
        }
    }
    return false;
}

using RmsdResidueKey = std::tuple<int, std::string, int, int>;

static std::map<RmsdResidueKey, std::vector<AtomCoordinate>> groupResiduesForRmsd(
    const std::vector<AtomCoordinate>& atoms
) {
    std::map<RmsdResidueKey, std::vector<AtomCoordinate>> residues;
    if (atoms.empty()) {
        return residues;
    }
    std::map<std::tuple<int, std::string, int, char, std::string>, int> occurrences;
    size_t residueStart = 0;
    for (size_t i = 1; i <= atoms.size(); i++) {
        bool endOfResidue = (i == atoms.size()) ||
                            startsNewResidue(atoms[i], atoms[i - 1]);
        if (!endOfResidue) {
            continue;
        }
        const AtomCoordinate& first = atoms[residueStart];
        std::tuple<int, std::string, int, char, std::string> baseKey(
            first.model, first.chain, first.residue_index, first.insertion_code, first.residue
        );
        RmsdResidueKey key(
            first.model, first.chain, first.residue_index, occurrences[baseKey]++
        );
        std::vector<AtomCoordinate>& residue = residues[key];
        residue.insert(residue.end(), atoms.begin() + residueStart, atoms.begin() + i);
        residueStart = i;
    }
    return residues;
}

static bool calculateAlignedRmsd(
    const std::vector<AtomCoordinate>& atoms1,
    const std::vector<AtomCoordinate>& atoms2,
    float& backboneRmsd,
    float& allAtomRmsd
) {
    backboneRmsd = 0.0f;
    allAtomRmsd = 0.0f;
    if (atoms1.empty() || atoms2.empty() || atoms1.size() != atoms2.size()) {
        return false;
    }

    std::vector<AtomCoordinate> alignedAtoms1;
    std::vector<AtomCoordinate> alignedAtoms2;
    std::vector<AtomCoordinate> alignedBackbone1;
    std::vector<AtomCoordinate> alignedBackbone2;

    std::map<RmsdResidueKey, std::vector<AtomCoordinate>> residues1 = groupResiduesForRmsd(atoms1);
    std::map<RmsdResidueKey, std::vector<AtomCoordinate>> residues2 = groupResiduesForRmsd(atoms2);
    if (residues1.size() != residues2.size()) {
        return false;
    }

    auto atomSorter = [](const AtomCoordinate& a, const AtomCoordinate& b) {
        if (a.atom != b.atom) {
            return a.atom < b.atom;
        }
        if (a.residue != b.residue) {
            return a.residue < b.residue;
        }
        if (a.chain != b.chain) {
            return a.chain < b.chain;
        }
        return a.atom_index < b.atom_index;
    };

    for (const auto& kv : residues1) {
        auto it = residues2.find(kv.first);
        if (it == residues2.end() || kv.second.empty() || it->second.empty()) {
            return false;
        }
        const AtomCoordinate& left = kv.second.front();
        const AtomCoordinate& right = it->second.front();
        if (left.residue != right.residue || kv.second.size() != it->second.size()) {
            return false;
        }

        std::vector<AtomCoordinate> sorted1 = kv.second;
        std::vector<AtomCoordinate> sorted2 = it->second;
        std::sort(sorted1.begin(), sorted1.end(), atomSorter);
        std::sort(sorted2.begin(), sorted2.end(), atomSorter);
        for (size_t atomIdx = 0; atomIdx < sorted1.size(); atomIdx++) {
            if (sorted1[atomIdx].atom != sorted2[atomIdx].atom) {
                return false;
            }
            alignedAtoms1.push_back(sorted1[atomIdx]);
            alignedAtoms2.push_back(sorted2[atomIdx]);
        }

        auto appendBackboneAtom = [&](const std::string& atomName) {
            auto leftIt = std::find_if(
                sorted1.begin(), sorted1.end(),
                [&](const AtomCoordinate& atom) { return atom.atom == atomName; }
            );
            auto rightIt = std::find_if(
                sorted2.begin(), sorted2.end(),
                [&](const AtomCoordinate& atom) { return atom.atom == atomName; }
            );
            if (leftIt != sorted1.end() && rightIt != sorted2.end()) {
                alignedBackbone1.push_back(*leftIt);
                alignedBackbone2.push_back(*rightIt);
            }
        };
        appendBackboneAtom("N");
        appendBackboneAtom("CA");
        appendBackboneAtom("C");
    }

    if (alignedAtoms1.empty() || alignedBackbone1.empty()) {
        return false;
    }
    backboneRmsd = RMSD(alignedBackbone1, alignedBackbone2);
    allAtomRmsd = RMSD(alignedAtoms1, alignedAtoms2);
    return true;
}

static std::string makeFragmentSuffix(const ContainerFragment& fragment, size_t index) {
    std::string suffix = "|m" + std::to_string(fragment.model);
    if (!fragment.chain.empty()) {
        suffix += "|c";
        suffix += fragment.chain;
    }
    suffix += "|f";
    suffix += std::to_string(index);
    return suffix;
}

int print_usage(void) {
    std::cout << "Usage: foldcomp compress <pdb|cif> [<fcz>]" << std::endl;
    std::cout << "       foldcomp compress [-t number] <dir|tar(.gz)> [<dir|tar|db>]" << std::endl;
    std::cout << "       foldcomp decompress <fcz|tar> [<pdb>]" << std::endl;
    std::cout << "       foldcomp decompress [-t number] <dir|tar(.gz)|db> [<dir|tar>]" << std::endl;
    std::cout << "       foldcomp extract [--plddt|--amino-acid] <fcz> [<fasta>]" << std::endl;
    std::cout << "       foldcomp extract [--plddt|--amino-acid] [-t number] <dir|tar(.gz)|db> [<fasta_out>]" << std::endl;
    std::cout << "       foldcomp check <fcz>" << std::endl;
    std::cout << "       foldcomp check [-t number] <dir|tar(.gz)|db>" << std::endl;
    std::cout << "       foldcomp rmsd <pdb|cif> <pdb|cif>" << std::endl;
    std::cout << " -h, --help               print this help message" << std::endl;
    std::cout << " -v, --version            print version" << std::endl;
    std::cout << " -t, --threads            threads for (de)compression of folders/tar files [default=1]" << std::endl;
    std::cout << " -r, --recursive          recursively look for files in directory [default=0]" << std::endl;
    std::cout << " -f, --file               input is a list of files [default=0]" << std::endl;
    std::cout << " -a, --alt                use alternative atom order [default=false]" << std::endl;
    std::cout << " -b, --break              interval size to save absolute atom coordinates [default=" << anchor_residue_threshold << "]" << std::endl;
    std::cout << " --max-backbone-rmsd      if set, spill FCZ fragments above this backbone RMSD to raw [default=inf]" << std::endl;
    std::cout << " -z, --tar                save as tar file [default=false]" << std::endl;
    std::cout << " -d, --db                 save as database [default=false]" << std::endl;
    std::cout << " -y, --overwrite          overwrite existing files [default=false]" << std::endl;
    std::cout << " -l, --id-list            a file of id list to be processed (only for database input)" << std::endl;
    std::cout << " -m, --id-mode            id mode for database input. 0: database keys, 1: names (.lookup) [default=1]" << std::endl;
    std::cout << " --skip-discontinuous     skip PDB with with discontinuous residues (only batch compression)" << std::endl;
    std::cout << " --check                  check FCZ before and skip entries with error (only for batch decompression)" << std::endl;
    std::cout << " --plddt                  extract pLDDT score (only for extraction mode)" << std::endl;
    std::cout << " -p, --plddt-digits       extract pLDDT score with specified number of digits (only for extraction mode)" << std::endl;
    std::cout << "                          - 1: single digit (fasta-like format), 2: 2-digit(00-99; tsv), 3: 3-digit, 4: 4-digit (max)" << std::endl;
    std::cout << " --fasta, --amino-acid    extract amino acid sequence (only for extraction mode)" << std::endl;
    std::cout << " --no-merge               do not merge output files (only for extraction mode)" << std::endl;
    std::cout << " --use-title              use TITLE as the output file name (only for extraction mode)" << std::endl;
    std::cout << " --time                   measure time for compression/decompression" << std::endl;
    std::cout << " --use-cache              use cached index for database input [default=false]" << std::endl;
    return 0;
}

inline int print_version(void) {
    std::cout << "foldcomp " << FOLDCOMP_VERSION << std::endl;
    return 0;
}

int rmsd(const std::string& pdb1, const std::string& pdb2) {
    // RMSD calculation between two PDB/mmCIF files
    StructureReader reader;
    reader.load(pdb1);
    std::vector<AtomCoordinate> atomCoordinates1;
    reader.readAllAtoms(atomCoordinates1);
    reader.load(pdb2);
    std::vector<AtomCoordinate> atomCoordinates2;
    reader.readAllAtoms(atomCoordinates2);
    removeAlternativePosition(atomCoordinates1);
    removeAlternativePosition(atomCoordinates2);
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

    float backboneRmsd = 0.0f;
    float allAtomRmsd = 0.0f;
    if (!calculateAlignedRmsd(atomCoordinates1, atomCoordinates2, backboneRmsd, allAtomRmsd)) {
        std::cerr << "[Error] Residue or atom identity mismatch between the two files." << std::endl;
        return 1;
    }
    // Print
    std::cout << pdb1 << '\t' << pdb2 << '\t';
    std::map<RmsdResidueKey, std::vector<AtomCoordinate>> residues1 = groupResiduesForRmsd(atomCoordinates1);
    std::cout << residues1.size() << '\t' << atomCoordinates1.size() << '\t';
    std::cout << backboneRmsd << '\t';
    std::cout << allAtomRmsd << std::endl;
    return 0;
}

int main(int argc, char* const *argv) {
    if (argc < 3) {
        // Check if version is requested
        if (argc == 2 && (strcmp(argv[1], "--version") == 0 || strcmp(argv[1], "-v") == 0)) {
            return print_version();
        }
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
    int check_before_decompression = 0;
    int id_mode = 1;
    int use_cache = 0;
    std::string user_id_file = "";
    std::vector<std::string> user_names;
    std::vector<uint32_t> user_ids;

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
            {"help",               no_argument,                           0, 'h'},
            {"alt",                no_argument,                           0, 'a'},
            {"tar",                no_argument,                           0, 'z'},
            {"recursive",          no_argument,                           0, 'r'},
            {"file",               no_argument,                           0, 'f'},
            {"plddt",              no_argument,                   &ext_mode,  0 },
            {"fasta",              no_argument,                   &ext_mode,  1 },
            {"amino-acid",         no_argument,                   &ext_mode,  1 },
            {"no-merge",           no_argument,                  &ext_merge,  0 },
            {"overwrite",          no_argument,                  &overwrite,  1 },
            {"time",               no_argument,               &measure_time,  1 },
            {"skip-discontinuous", no_argument,         &skip_discontinuous,  1 },
            {"check",              no_argument, &check_before_decompression,  1 },
            {"use-title",          no_argument,              &ext_use_title,  1 },
            {"db",                 no_argument,                           0, 'd'},
            {"version",            no_argument,                           0, 'v'},
            {"threads",      required_argument,                           0, 't'},
            {"break",        required_argument,                           0, 'b'},
            {"max-backbone-rmsd", required_argument,                      0, 1000},
            {"id-list",      required_argument,                           0, 'l'},
            {"id-mode",      required_argument,                           0, 'm'},
            {"plddt-digits", required_argument,                           0, 'p'},
            {"use-cache",      no_argument,                          &use_cache,  1 },
            {0,                              0,                           0,  0 }
    };

    // Parse command line options with getopt_long
    int flag = getopt_long(argc, argv, "hadzrfyvt:b:l:p:", long_options, &option_index);
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
            case 'y':
                overwrite = 1;
                break;
            case 'b':
                anchor_residue_threshold = atoi(optarg);
                break;
            case 1000:
                max_backbone_rmsd = static_cast<float>(atof(optarg));
                break;
            case 'l':
                user_id_file = std::string(optarg);
                break;
            case 'm':
                id_mode = atoi(optarg);
                if (id_mode != 0 && id_mode != 1) {
                    std::cerr << "[Error] Invalid id mode. Please use 0 or 1." << std::endl;
                    return print_usage();
                }
                break;
            case 'd':
                db_output = 1;
                break;
            case 'p':
                ext_plddt_digits = atoi(optarg);
                break;
            case 'v':
                return print_version();
            case '?':
                return print_usage();
            default:
                break;
        }
        flag = getopt_long(argc, argv, "hadzrfyt:b:l:p:", long_options, &option_index);
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
            outputSuffix = "plddt";
            if (ext_plddt_digits != 1) {
                outputSuffix = "plddt.tsv";
            }
        } else if (ext_mode == 1) {
            outputSuffix = "fasta";
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
    std::vector<std::string> single_file_inputs;
    if (file_input) {
        FILE* inputFile = fopen(input.c_str(), "rb");
        if (inputFile == nullptr) {
            std::cerr << "[Error] Could not open file " << input << std::endl;
            return EXIT_FAILURE;
        }
        char lineBuffer[4096];
        while (fgets(lineBuffer, sizeof(lineBuffer), inputFile) != nullptr) {
            std::string line(lineBuffer);
            while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
                line.pop_back();
            }
            // If the file ends with .pdb, .pdb.gz, .cif, .cif.gz, .fcz, assume it is a single file input
            if (stringEndsWith(".pdb", line) || stringEndsWith(".pdb.gz", line) ||
                stringEndsWith(".cif", line) || stringEndsWith(".cif.gz", line) ||
                stringEndsWith(".fcz", line)) {
                single_file_inputs.push_back(line);
            } else {
                inputs.push_back(line);
            }
        }
        fclose(inputFile);
    } else {
        inputs.push_back(input);
    }

    std::string output;
    bool useMMCIFOutput = false;
    if (argc == optind + 3) {
        has_output = 1;
        output = argv[optind + 2];
        if (stringEndsWith(".tar", output)) {
            save_as_tar = 1;
        }
        if (mode == DECOMPRESS && isMMCIFOutputPath(output)) {
            useMMCIFOutput = true;
            outputSuffix = stringEndsWith(".mmcif", output) || stringEndsWith(".mmcif.tar", output) ? "mmcif" : "cif";
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

        unsigned int key = 0;
        for (size_t i = 0; i < inputs.size() + 1; i++) {
            const std::string& input = (i == inputs.size()) ? "" : inputs[i];
            Processor* processor;
            struct stat st;
            if (i != inputs.size()) {
                if (stringEndsWith(".tar", input) || stringEndsWith(".tar.gz", input) || stringEndsWith(".tgz", input)) {
                    processor = new TarProcessor(input);
                }
                else if (stat((input + ".dbtype").c_str(), &st) == 0) {
                    if (user_id_file.size() > 0) {
                        processor = new DatabaseProcessor(input, user_id_file, id_mode, use_cache);
                    }
                    else {
                        processor = new DatabaseProcessor(input);
                    }
                }
#ifdef HAVE_GCS
                else if (stringStartsWith("gcs://", input)) {
                    processor = new GcsProcessor(input);
                }
#endif
                else {
                    processor = new DirectoryProcessor(input, recursive);
                }
            } else {
                if (single_file_inputs.size() > 0) {
                    processor = new DirectoryProcessor(single_file_inputs);
                } else {
                    continue;
                }
            }
            process_entry_func func = [&](const char* name, const char* dataBuffer, size_t size) -> bool {
                TimerGuard guard(name, measure_time);
                std::vector<AtomCoordinate> atomCoordinates;
                std::vector<BackboneChain> compData;
                StructureReader reader;

                std::string base = baseName(name);
                std::pair<std::string, std::string> outputParts = getFileParts(base);
                std::string outputFile;
                if (save_as_tar) {
                    outputFile = outputParts.first;
                } else if (db_output) {
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
                    std::cerr << "[Warning] No protein atoms found in the input file: " << base << std::endl;
                    return true;
                }

                // replace the title with only the base name if it was the same as the file name
                std::string title = reader.title == base ? outputParts.first : reader.title;

                std::vector<ContainerFragment> encodedFragments;
                std::vector<std::pair<size_t, size_t>> chain_indices = identifyChains(atomCoordinates);
                // Check if there are multiple chains or regions with discontinous residue indices
                for (size_t i = 0; i < chain_indices.size(); i++) {
                    std::vector<std::pair<size_t, size_t>> frag_indices = identifyDiscontinousResInd(
                        atomCoordinates, chain_indices[i].first, chain_indices[i].second
                    );
                    if (frag_indices.empty()) {
                        frag_indices.push_back(chain_indices[i]);
                    }
                    if (skip_discontinuous && frag_indices.size() > 1) {
                        std::string message = "Skipping discontinuous chain: " + base + "\n";
                        std::cerr << message;
                        continue;
                    }
                    for (size_t j = 0; j < frag_indices.size(); j++) {
                        tcb::span<AtomCoordinate> chain_span = tcb::span<AtomCoordinate>(
                            &atomCoordinates[frag_indices[j].first],
                            atomCoordinates.data() + frag_indices[j].second
                        );
                        std::vector<BackboneRegion> regions = identifyBackboneRegions(chain_span);
                        if (regions.empty()) {
                            std::cerr << "[Warning] Skipping empty fragment: " << base << std::endl;
                            continue;
                        }
                        std::vector<bool> region_needs_raw(regions.size(), false);
                        for (size_t region_index = 0; region_index < regions.size(); region_index++) {
                            if (!regions[region_index].encodable) {
                                continue;
                            }
                            tcb::span<AtomCoordinate> region_span(
                                chain_span.data() + regions[region_index].start,
                                chain_span.data() + regions[region_index].end
                            );
                            region_needs_raw[region_index] = regionNeedsRawFallback(region_span);
                        }
                        if (shouldStoreSmallMixedFragmentAsRaw(chain_span, regions, region_needs_raw)) {
                            std::cerr << "[Warning] Using raw encoding for short mixed fragment: "
                                      << base << std::endl;
                            ContainerFragment containerFragment;
                            containerFragment.kind = CONTAINER_FRAGMENT_KIND_RAW_ATOMS;
                            containerFragment.model = chain_span.front().model;
                            containerFragment.chain = chain_span.front().chain;
                            if (!serializeAtomCoordinates(
                                    tcb::span<const AtomCoordinate>(chain_span.data(), chain_span.size()),
                                    containerFragment.payload)) {
                                std::cerr << "[Error] Failed to serialize raw fragment: " << base << std::endl;
                                return false;
                            }
                            encodedFragments.push_back(std::move(containerFragment));
                            continue;
                        }
                        for (size_t region_index = 0; region_index < regions.size(); region_index++) {
                            const auto& region = regions[region_index];
                            tcb::span<AtomCoordinate> region_span(
                                chain_span.data() + region.start,
                                chain_span.data() + region.end
                            );
                            ContainerFragment containerFragment;
                            containerFragment.model = chain_span[region.start].model;
                            containerFragment.chain = chain_span[region.start].chain;
                            bool useFoldcompEncoding = region.encodable && !region_needs_raw[region_index];
                            if (useFoldcompEncoding) {
                                Foldcomp compRes;
                                compRes.strTitle = title;
                                compRes.anchorThreshold = anchor_residue_threshold;
                                compData = compRes.compress(region_span);
                                if (compData.empty()) {
                                    std::cerr << "[Warning] Skipping fragment with incomplete backbone: "
                                              << base << std::endl;
                                    continue;
                                }
                                containerFragment.kind = CONTAINER_FRAGMENT_KIND_FCZ;
                                if (compRes.writeString(containerFragment.payload) != 0) {
                                    std::cerr << "[Error] Failed to serialize fragment payload: " << base << std::endl;
                                    return false;
                                }
                                compData.clear();
                                if (compRes.exceedsBackboneRmsdThreshold(max_backbone_rmsd)) {
                                    std::cerr << "[Warning] Falling back to raw due to backbone RMSD threshold: "
                                              << base << std::endl;
                                    containerFragment.kind = CONTAINER_FRAGMENT_KIND_RAW_ATOMS;
                                    containerFragment.payload.clear();
                                    if (!serializeAtomCoordinates(
                                            tcb::span<const AtomCoordinate>(region_span.data(), region_span.size()),
                                            containerFragment.payload)) {
                                        std::cerr << "[Error] Failed to serialize raw fragment: " << base << std::endl;
                                        return false;
                                    }
                                }
                            } else {
                                containerFragment.kind = CONTAINER_FRAGMENT_KIND_RAW_ATOMS;
                                if (!serializeAtomCoordinates(
                                        tcb::span<const AtomCoordinate>(region_span.data(), region_span.size()),
                                        containerFragment.payload)) {
                                    std::cerr << "[Error] Failed to serialize raw fragment: " << base << std::endl;
                                    return false;
                                }
                            }
                            encodedFragments.push_back(std::move(containerFragment));
                        }
                    }
                }

                if (encodedFragments.empty()) {
                    std::cerr << "[Warning] No protein fragments found in input file: " << base << std::endl;
                    return true;
                }

                bool useContainer = encodedFragments.size() > 1 ||
                                    encodedFragments[0].kind != CONTAINER_FRAGMENT_KIND_FCZ ||
                                    encodedFragments[0].chain.size() != 1 ||
                                    encodedFragments[0].model != 1;
                std::string encoded;
                if (useContainer) {
                    if (!writeContainerToString(encoded, title, encodedFragments)) {
                        std::cerr << "[Error] Failed to write container payload: " << base << std::endl;
                        return false;
                    }
                } else {
                    encoded = encodedFragments[0].payload;
                }

                std::string filename = outputFile;
                if (!db_output) {
                    if (isCompressible(outputParts)) {
                        filename += ".fcz";
                    } else {
                        filename += "." + outputParts.second;
                    }
                }

                if (db_output) {
                    std::string dbKey = baseName(filename);
                    std::replace(dbKey.begin(), dbKey.end(), '.', '_');
#pragma omp critical
                    {
                        writer_append(handle, encoded.c_str(), encoded.size(), key, dbKey.c_str());
                        key++;
                    }
                } else if (save_as_tar) {
#pragma omp critical
                    {
                        mtar_write_file_header(&tar_out, baseName(filename).c_str(), encoded.size());
                        mtar_write_data(&tar_out, encoded.c_str(), encoded.size());
                    }
                } else {
                    if (stat(filename.c_str(), &st) == 0 && !overwrite) {
                        std::cerr << "[Error] Output file already exists: " << baseName(outputFile) << std::endl;
                        return false;
                    }
                    FILE* outFile = fopen(filename.c_str(), "wb");
                    if (outFile == nullptr) {
                        std::cerr << "[Error] Failed to open output file: " << filename << std::endl;
                        return false;
                    }
                    size_t written = fwrite(encoded.data(), 1, encoded.size(), outFile);
                    fclose(outFile);
                    if (written != encoded.size()) {
                        std::cerr << "[Error] Failed to write output file: " << filename << std::endl;
                        return false;
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

        unsigned int key = 0;
        for (size_t i = 0; i < inputs.size() + 1; i++) {
            const std::string& input = (i == inputs.size()) ? "" : inputs[i];
            Processor* processor;
            struct stat st;
            if (i != inputs.size()) {
                if (stringEndsWith(".tar", input) || stringEndsWith(".tar.gz", input) || stringEndsWith(".tgz", input)) {
                    processor = new TarProcessor(input);
                }
                else if (stat((input + ".dbtype").c_str(), &st) == 0) {
                    if (user_id_file.size() > 0) {
                        processor = new DatabaseProcessor(input, user_id_file, id_mode, use_cache);
                    }
                    else {
                        processor = new DatabaseProcessor(input);
                    }
                }
#ifdef HAVE_GCS
                else if (stringStartsWith("gcs://", input)) {
                    processor = new GcsProcessor(input);
                }
#endif
                else {
                    processor = new DirectoryProcessor(input, recursive);
                }
            } else {
                if (single_file_inputs.size() > 0) {
                    processor = new DirectoryProcessor(single_file_inputs);
                }
                else {
                    continue;
                }
            }

            process_entry_func func = [&](const char* name, const char* dataBuffer, size_t size) -> bool {
                TimerGuard guard(name, measure_time);
                std::vector<DecompressedSegment> segments;
                std::string structureTitle;

                std::string containerTitle;
                std::vector<ContainerFragment> containerFragments;
                bool isContainer = readContainer(dataBuffer, size, containerTitle, containerFragments);
                if (isContainer) {
                    structureTitle = containerTitle;
                    for (const auto& fragment : containerFragments) {
                        std::vector<AtomCoordinate> atomCoordinates;
                        if (fragment.kind == CONTAINER_FRAGMENT_KIND_RAW_ATOMS) {
                            if (!deserializeAtomCoordinates(
                                    fragment.payload.data(), fragment.payload.size(), atomCoordinates)) {
                                std::cerr << "[Error] Failed to decode a raw atom fragment." << std::endl;
                                return false;
                            }
                        } else if (fragment.kind == CONTAINER_FRAGMENT_KIND_FCZ) {
                            Foldcomp compRes;
                            int flag = compRes.read(fragment.payload.data(), fragment.payload.size());
                            if (flag != 0) {
                                std::cerr << "[Error] Failed to read a container fragment." << std::endl;
                                return false;
                            }
                            compRes.useAltAtomOrder = use_alt_order;
                            if (check_before_decompression) {
                                ValidityError err = compRes.checkValidity();
                                if (err != ValidityError::SUCCESS) {
                                    printValidityError(err, compRes.strTitle);
                                    return true;
                                }
                            }
                            flag = compRes.decompress(atomCoordinates);
                            if (flag != 0) {
                                std::cerr << "[Error] decompressing container fragment." << std::endl;
                                return false;
                            }
                            if (structureTitle.empty()) {
                                structureTitle = compRes.strTitle;
                            }
                        } else {
                            std::cerr << "[Error] Unsupported container fragment kind: "
                                      << static_cast<int>(fragment.kind) << std::endl;
                            return false;
                        }
                        for (auto& atom : atomCoordinates) {
                            atom.model = fragment.model;
                            atom.chain = fragment.chain;
                        }
                        segments.push_back({fragment.model, std::move(atomCoordinates)});
                    }
                } else {
                    Foldcomp compRes;
                    int flag = compRes.read(dataBuffer, size);
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
                    compRes.useAltAtomOrder = use_alt_order;
                    if (check_before_decompression) {
                        ValidityError err = compRes.checkValidity();
                        if (err != ValidityError::SUCCESS) {
                            printValidityError(err, compRes.strTitle);
                            return true;
                        }
                    }
                    std::vector<AtomCoordinate> atomCoordinates;
                    flag = compRes.decompress(atomCoordinates);
                    if (flag != 0) {
                        std::cerr << "[Error] decompressing compressed data." << std::endl;
                        return false;
                    }
                    segments.push_back({1, std::move(atomCoordinates)});
                    structureTitle = compRes.strTitle;
                }

                if (segments.empty()) {
                    std::cerr << "[Error] No segments were decompressed." << std::endl;
                    return false;
                }

                std::string base = baseName(name);
                std::pair<std::string, std::string> outputParts = getFileParts(base);
                std::string outputFile;
                if (save_as_tar) {
                    outputFile = outputParts.first + "." + outputSuffix;
                } else if (db_output) {
                    outputFile = outputParts.first;
                } else if (isSingleFileInput) {
                    outputFile = output;
                } else {
                    outputFile = output + "/" + outputParts.first + "." + outputSuffix;
                }

                std::string structureText;
#ifdef FOLDCOMP_WITH_MMCIF_OUTPUT
                if (useMMCIFOutput) {
                    writeSegmentsToMMCIF(segments, structureTitle, structureText);
                } else {
                    writeSegmentsToPDB(segments, structureTitle, structureText);
                }
#else
                writeSegmentsToPDB(segments, structureTitle, structureText);
#endif

                if (db_output) {
                    structureText.push_back('\0');
#pragma omp critical
                    {
                        writer_append(handle, structureText.c_str(), structureText.size(), key, outputFile.c_str());
                        key++;
                    }
                } else if (save_as_tar) {
#pragma omp critical
                    {
                        mtar_write_file_header(&tar_out, outputFile.c_str(), structureText.size());
                        mtar_write_data(&tar_out, structureText.c_str(), structureText.size());
                    }
                } else {
                    // Write decompressed data to file
                    // Check output file exists
                    if (stat(outputFile.c_str(), &st) == 0 && !overwrite) {
                        std::cerr << "[Error] Output file already exists: " << baseName(outputFile) << std::endl;
                        return false;
                    }
                    FILE* out = fopen(outputFile.c_str(), "wb");
                    if (out == nullptr) {
                        std::cerr << "[Error] Writing decompressed data to file: " << output << std::endl;
                        return false;
                    }
                    size_t written = fwrite(structureText.data(), 1, structureText.size(), out);
                    fclose(out);
                    if (written != structureText.size()) {
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
            if (stat(output.c_str(), &st) == -1 && !ext_merge) {
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
            } else if (!ext_merge) {
                std::cout << "Output directory: " << output << std::endl;
            } else {
                // Single file output. Remove "/" from end of output
                if (output[output.size() - 1] == '/') {
                    output = output.substr(0, output.size() - 1);
                }
                std::cout << "Output: " << output << std::endl;
            }
        }

        bool isMergedOutput = false;
        FILE* default_out = nullptr;
        if (!save_as_tar && !db_output && !isSingleFileInput && ext_merge) {
            default_out = fopen(output.c_str(), "wb");
            isMergedOutput = true;
        }

        unsigned int key = 0;
        for (size_t i = 0; i < inputs.size() + 1; i++) {
            const std::string& input = (i == inputs.size()) ? "" : inputs[i];
            Processor* processor;
            struct stat st;
            if (i != inputs.size()) {
                if (stringEndsWith(".tar", input) || stringEndsWith(".tar.gz", input) || stringEndsWith(".tgz", input)) {
                    processor = new TarProcessor(input);
                }
                else if (stat((input + ".dbtype").c_str(), &st) == 0) {
                    if (user_id_file.size() > 0) {
                        processor = new DatabaseProcessor(input, user_id_file, id_mode, use_cache);
                    }
                    else {
                        processor = new DatabaseProcessor(input);
                    }
                }
#ifdef HAVE_GCS
                else if (stringStartsWith("gcs://", input)) {
                    processor = new GcsProcessor(input);
                }
#endif
                else {
                    processor = new DirectoryProcessor(input, recursive);
                }
            } else {
                if (single_file_inputs.size() > 0) {
                    processor = new DirectoryProcessor(single_file_inputs);
                }
                else {
                    continue;
                }
            }

            process_entry_func func = [&](const char* name, const char* dataBuffer, size_t size) -> bool {
                std::string strName(name);
                std::vector<ContainerFragment> containerFragments;
                std::string containerTitle;
                bool isContainer = readContainer(dataBuffer, size, containerTitle, containerFragments);

                std::string extractedText;
                if (!isContainer) {
                    std::string outputTitle = strName;
                    std::string extractedData;
                    int totalResidue = 0;
                    Foldcomp compRes;
                    int flag = compRes.read(dataBuffer, size);
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
                    if (ext_use_title && !compRes.strTitle.empty()) {
                        outputTitle = compRes.strTitle;
                    }
                    compRes.extract(extractedData, ext_mode, ext_plddt_digits);
                    totalResidue = compRes.nResidue;

                    if (ext_mode == 0 && ext_plddt_digits > 1) {
                        extractedText += outputTitle + "\t" + std::to_string(totalResidue) + "\t" + extractedData + "\n";
                    } else {
                        extractedText += ">" + outputTitle + "\n" + extractedData + "\n";
                    }
                } else {
                    struct ExtractGroup {
                        int16_t model = 1;
                        std::string chain;
                        int totalResidue = 0;
                        std::string data;
                    };
                    std::vector<ExtractGroup> groups;
                    std::string baseTitle = ext_use_title ? containerTitle : strName;

                    for (size_t fragmentIndex = 0; fragmentIndex < containerFragments.size(); fragmentIndex++) {
                        Foldcomp compRes;
                        int flag = compRes.read(containerFragments[fragmentIndex].payload.data(),
                                                containerFragments[fragmentIndex].payload.size());
                        if (flag != 0) {
                            std::cerr << "[Error] Failed to read a container fragment during extraction." << std::endl;
                            return false;
                        }
                        if (baseTitle.empty() && ext_use_title && !compRes.strTitle.empty()) {
                            baseTitle = compRes.strTitle;
                        }

                        size_t groupIdx = groups.size();
                        for (size_t i = 0; i < groups.size(); i++) {
                            if (groups[i].model == containerFragments[fragmentIndex].model &&
                                groups[i].chain == containerFragments[fragmentIndex].chain) {
                                groupIdx = i;
                                break;
                            }
                        }
                        if (groupIdx == groups.size()) {
                            ExtractGroup g;
                            g.model = containerFragments[fragmentIndex].model;
                            g.chain = containerFragments[fragmentIndex].chain;
                            groups.push_back(std::move(g));
                        }

                        std::string fragmentData;
                        compRes.extract(fragmentData, ext_mode, ext_plddt_digits);
                        if (ext_mode == 0 && ext_plddt_digits > 1 && !groups[groupIdx].data.empty()) {
                            groups[groupIdx].data.push_back(',');
                        }
                        groups[groupIdx].data += fragmentData;
                        groups[groupIdx].totalResidue += compRes.nResidue;
                    }

                    if (baseTitle.empty()) {
                        baseTitle = strName;
                    }
                    bool annotateGroup = groups.size() > 1;
                    for (const auto& group : groups) {
                        std::string groupTitle = baseTitle;
                        if (annotateGroup) {
                            std::string suffix = "|m" + std::to_string(group.model);
                            if (!group.chain.empty()) {
                                suffix += "|c";
                                suffix += group.chain;
                            }
                            groupTitle += suffix;
                        }
                        if (ext_mode == 0 && ext_plddt_digits > 1) {
                            extractedText += groupTitle + "\t" + std::to_string(group.totalResidue) + "\t" + group.data + "\n";
                        } else {
                            extractedText += ">" + groupTitle + "\n" + group.data + "\n";
                        }
                    }
                }
                std::string base = baseName(name);
                std::pair<std::string, std::string> outputParts = getFileParts(base);
                std::string outputFile;
                if (save_as_tar) {
                    outputFile = outputParts.first + "." + outputSuffix;
                } else if (db_output) {
                    outputFile = outputParts.first;
                } else if (isSingleFileInput) {
                    outputFile = output;
                } else {
                    outputFile = output + "/" + outputParts.first + "." + outputSuffix;
                }

                if (isMergedOutput) {
#pragma omp critical
                    {
                        fwrite(extractedText.data(), 1, extractedText.size(), default_out);
                    }
                } else if (db_output) {
                    extractedText.push_back('\0');
#pragma omp critical
                    {
                        writer_append(handle, extractedText.c_str(), extractedText.size(), key, outputFile.c_str());
                        key++;
                    }
                } else if (save_as_tar) {
#pragma omp critical
                    {
                        mtar_write_file_header(&tar_out, outputFile.c_str(), extractedText.size());
                        mtar_write_data(&tar_out, extractedText.c_str(), extractedText.size());
                    }
                } else {
                    FILE* outputFileHandle = fopen(outputFile.c_str(), "wb");
                    if (outputFileHandle == nullptr) {
                        std::cerr << "[Error] Could not open file " << outputFile << std::endl;
                        return false;
                    }
                    size_t written = fwrite(extractedText.data(), 1, extractedText.size(), outputFileHandle);
                    fclose(outputFileHandle);
                    if (written != extractedText.size()) {
                        std::cerr << "[Error] Could not write file " << outputFile << std::endl;
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
        } else if (isMergedOutput && default_out != nullptr) {
            fclose(default_out);
        }
    } else if (mode == CHECK) {
        if (inputs.size() == 1) {
            std::cout << "Checking " << input << std::endl;
        } else {
            std::cout << "Checking files in " << input;
            std::cout << " using " << num_threads << " threads" << std::endl;
        }

        for (size_t i = 0; i < inputs.size() + 1; i++) {
            const std::string& input = (i == inputs.size()) ? "" : inputs[i];
            Processor* processor;
            struct stat st;
            if (i != inputs.size()) {
                if (stringEndsWith(".tar", input) || stringEndsWith(".tar.gz", input) || stringEndsWith(".tgz", input)) {
                    processor = new TarProcessor(input);
                }
                else if (stat((input + ".dbtype").c_str(), &st) == 0) {
                    if (user_id_file.size() > 0) {
                        processor = new DatabaseProcessor(input, user_id_file, id_mode, use_cache);
                    }
                    else {
                        processor = new DatabaseProcessor(input);
                    }
                }
#ifdef HAVE_GCS
                else if (stringStartsWith("gcs://", input)) {
                    processor = new GcsProcessor(input);
                }
#endif
                else {
                    processor = new DirectoryProcessor(input, recursive);
                }
            } else {
                if (single_file_inputs.size() > 0) {
                    processor = new DirectoryProcessor(single_file_inputs);
                }
                else {
                    continue;
                }
            }

            process_entry_func func = [&](const char* name, const char* dataBuffer, size_t size) -> bool {
                std::vector<ContainerFragment> containerFragments;
                std::string containerTitle;
                bool isContainer = readContainer(dataBuffer, size, containerTitle, containerFragments);
                if (isContainer) {
                    for (size_t fragmentIndex = 0; fragmentIndex < containerFragments.size(); fragmentIndex++) {
                        if (containerFragments[fragmentIndex].kind == CONTAINER_FRAGMENT_KIND_RAW_ATOMS) {
                            continue;
                        }
                        if (containerFragments[fragmentIndex].kind != CONTAINER_FRAGMENT_KIND_FCZ) {
                            std::cerr << "[Error] Unsupported container fragment kind during check: "
                                      << static_cast<int>(containerFragments[fragmentIndex].kind) << std::endl;
                            return false;
                        }
                        Foldcomp compRes;
                        int flag = compRes.read(containerFragments[fragmentIndex].payload.data(),
                                                containerFragments[fragmentIndex].payload.size());
                        if (flag != 0) {
                            std::cerr << "[Error] Failed to read a container fragment during check." << std::endl;
                            return false;
                        }
                        ValidityError err = compRes.checkValidity();
                        std::string sname(name);
                        sname += makeFragmentSuffix(containerFragments[fragmentIndex], fragmentIndex);
                        printValidityError(err, sname);
                    }
                    return true;
                } else {
                    Foldcomp compRes;
                    int flag = compRes.read(dataBuffer, size);
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
                }
            };
            processor->run(func, num_threads);
            delete processor;
        }
    }
    return EXIT_SUCCESS;
}
