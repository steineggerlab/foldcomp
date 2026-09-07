#pragma once

#include "atom_coordinate.h"

#include <cstdint>
#include <string>
#include <vector>

static constexpr size_t CONTAINER_MAGIC_LENGTH = 4;
static constexpr char CONTAINER_MAGIC[] = "FCZC";
static constexpr uint8_t FOLDCOMP_FORMAT_VERSION_CONTAINER = 2;
static constexpr uint8_t FOLDCOMP_FORMAT_FLAG_CONTAINER = 1;
static constexpr uint8_t CONTAINER_FRAGMENT_KIND_FCZ = 0;
static constexpr uint8_t CONTAINER_FRAGMENT_KIND_RAW_ATOMS = 2;

struct ContainerFragment {
    uint8_t kind = CONTAINER_FRAGMENT_KIND_FCZ;
    int16_t model = 1;
    std::string chain;
    std::string payload;
};

struct DecompressedSegment {
    int model = 1;
    std::vector<AtomCoordinate> atoms;
};

enum ParsePdbStatus {
    PARSE_PDB_OK = 0,
    PARSE_PDB_NO_ATOM = 1,
    PARSE_PDB_MULTIPLE_CHAINS = 2,
    PARSE_PDB_INVALID_FORMAT = 3,
};

bool hasContainerMagic(const char* data, size_t size);
bool parseStructureAtoms(
    const char* data, size_t size, bool singleChainOnly,
    std::vector<AtomCoordinate>& atoms, int& status,
    std::string* title = nullptr, const char* format = nullptr
);
int encodeStructureToFoldcomp(
    const std::string& title,
    const char* data,
    size_t size,
    const char* format,
    int anchorResidueThreshold,
    float maxBackboneRmsd,
    std::string& output
);
bool parsePDBAtoms(
    const char* data, size_t size, bool singleChainOnly,
    std::vector<AtomCoordinate>& atoms, int& status
);
bool readContainer(
    const char* data, size_t size, std::string& title,
    std::vector<ContainerFragment>& fragments
);
bool writeContainerToString(
    std::string& output, const std::string& title,
    const std::vector<ContainerFragment>& fragments
);
bool decodeStructureToSegments(
    const char* data, size_t size, bool useAltOrder, std::string& title,
    std::vector<DecompressedSegment>& segments
);
bool decodeStructureToAtoms(
    const char* data, size_t size, bool useAltOrder, std::string& title,
    std::vector<AtomCoordinate>& atoms
);
bool decodeStructureToPDB(
    const char* data, size_t size, bool useAltOrder, std::string& title,
    std::string& pdbText
);
// `precomputedLines` (optional): PDB_ATOM_LINE_LEN-byte ATOM lines for every atom
// across all segments, in segment/atom order (as produced by the GPU PDB writer).
// When provided, ATOM records are copied from it instead of formatted on the host.
void writeSegmentsToPDB(
    const std::vector<DecompressedSegment>& segments,
    const std::string& title,
    std::string& output,
    const char* precomputedLines = nullptr
);

// Per-segment metadata for the fully line-driven PDB assembler. A segment's model
// and chain are uniform across its atoms (as the GPU pipeline produces them), so
// only these three fields are needed to place MODEL/ENDMDL/TER records; all ATOM
// record content (including the TER fields) comes from `lines`. `chain` holds the
// full ChainId (not just the PDB column-21 byte) so segment-boundary comparisons
// match writeSegmentsToPDB's full-width chain comparison.
struct PdbLineSegment {
    int model;
    ChainId chain;
    int n_atoms;
};

// Assemble a byte-exact PDB from GPU-formatted ATOM lines plus per-segment
// metadata — no per-atom AtomCoordinate objects. `lines` holds
// sum(n_atoms) * PDB_ATOM_LINE_LEN bytes in segment/atom order. Produces output
// identical to writeSegmentsToPDB.
void writeLineSegmentsToPDB(
    const std::vector<PdbLineSegment>& segments,
    const std::string& title,
    const char* lines,
    std::string& output
);

#ifdef FOLDCOMP_WITH_MMCIF_OUTPUT
bool decodeStructureToMMCIF(
    const char* data, size_t size, bool useAltOrder, std::string& title,
    std::string& mmcifText
);
void writeSegmentsToMMCIF(
    const std::vector<DecompressedSegment>& segments,
    const std::string& title,
    std::string& output
);
#endif
