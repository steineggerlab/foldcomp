/**
 * File: atom_coordinate.h
 * Project: foldcomp
 * Created: 2021-01-18 12:43:08
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     The data type to handle atom coordinate comes here.
 * ---
 * Last Modified: 2022-11-29 14:39:06
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
#pragma once
#include "float3d.h"
#include "tcbspan.h"

#include <cstdint>
#include <cstring>
#include <fstream>
#include <ostream>
#include <string>
#include <vector>

/**
 * Fixed-length string stored as char array — zero heap allocation, minimal copy cost.
 * N is the maximum number of characters (excluding null terminator).
 * Supports implicit conversion to std::string and comparison with const char* / std::string.
 */
template<int N>
struct FixedStr {
    char data[N + 1];
    FixedStr() { data[0] = '\0'; }
    FixedStr(const char* s) { strncpy(data, s, N); data[N] = '\0'; }
    FixedStr(const std::string& s) { strncpy(data, s.c_str(), N); data[N] = '\0'; }
    FixedStr& operator=(const char* s)         { strncpy(data, s, N); data[N] = '\0'; return *this; }
    FixedStr& operator=(const std::string& s)  { strncpy(data, s.c_str(), N); data[N] = '\0'; return *this; }
    bool operator==(const char* s)             const { return strcmp(data, s) == 0; }
    bool operator==(const std::string& s)      const { return strcmp(data, s.c_str()) == 0; }
    bool operator==(const FixedStr& o)         const { return strcmp(data, o.data) == 0; }
    bool operator!=(const char* s)             const { return strcmp(data, s) != 0; }
    bool operator!=(const std::string& s)      const { return strcmp(data, s.c_str()) != 0; }
    bool operator!=(const FixedStr& o)         const { return strcmp(data, o.data) != 0; }
    bool operator< (const FixedStr& o)         const { return strcmp(data, o.data) <  0; }
    char& operator[](size_t i)       { return data[i]; }
    char  operator[](size_t i) const { return data[i]; }
    size_t size()      const { return strlen(data); }
    bool empty()       const { return data[0] == '\0'; }
    const char* c_str() const { return data; }
    operator std::string() const { return std::string(data); }
    template<int M> friend std::ostream& operator<<(std::ostream& os, const FixedStr<M>& s);
};
template<int N>
inline std::ostream& operator<<(std::ostream& os, const FixedStr<N>& s) { return os << s.data; }
// String concatenation: FixedStr + char* / string and vice versa
template<int N> inline std::string operator+(const FixedStr<N>& a, const char* b)         { return std::string(a.data) + b; }
template<int N> inline std::string operator+(const char* b, const FixedStr<N>& a)         { return b + std::string(a.data); }
template<int N> inline std::string operator+(const FixedStr<N>& a, const std::string& b)  { return std::string(a.data) + b; }
template<int N> inline std::string operator+(const std::string& b, const FixedStr<N>& a)  { return b + std::string(a.data); }
template<int N, int M> inline std::string operator+(const FixedStr<N>& a, const FixedStr<M>& b) { return std::string(a.data) + b.data; }
// Reverse comparison operators
template<int N> inline bool operator==(const char* s, const FixedStr<N>& a)        { return a == s; }
template<int N> inline bool operator!=(const char* s, const FixedStr<N>& a)        { return a != s; }
template<int N> inline bool operator==(const std::string& s, const FixedStr<N>& a) { return a == s; }
template<int N> inline bool operator!=(const std::string& s, const FixedStr<N>& a) { return a != s; }

using AtomName    = FixedStr<4>;  // e.g. "N", "CA", "OXT"  — max 4 chars
using ResidueName = FixedStr<3>;  // e.g. "ALA", "GLY"      — max 3 chars

// mmCIF auth_asym_id limit; also the width GPU-side chain buffers (chain_per_struct,
// PendingOXTState::chain) are packed to, so this must stay in sync with ChainId below.
constexpr int CHAIN_ID_LENGTH = 4;
using ChainId = FixedStr<CHAIN_ID_LENGTH>;  // e.g. "A", "AB"

class AtomCoordinate {
public:
    AtomCoordinate() = default;
    AtomCoordinate(
        AtomName a, ResidueName r, ChainId c,
        int ai, int ri, float x, float y, float z,
        float occupancy = 0.0f, float tempFactor = 0.0f,
        int model = 1, char insertionCode = ' ', char altloc = ' '
    );
    AtomCoordinate(
        AtomName a, ResidueName r, ChainId c,
        int ai, int ri, float3d coord,
        float occupancy = 0.0f, float tempFactor = 0.0f,
        int model = 1, char insertionCode = ' ', char altloc = ' '
    );
    // data
    AtomName    atom;
    ResidueName residue;
    ChainId     chain;
    int atom_index;
    int residue_index;
    float3d coordinate;
    float occupancy;
    float tempFactor;
    int model = 1;
    char insertion_code = ' ';
    char altloc = ' ';

    // operators
    bool operator==(const AtomCoordinate& other) const;
    bool operator!=(const AtomCoordinate& other) const;
    //BackboneChain toCompressedResidue();

    //methods
    bool isBackbone() const;
    void print(int option = 0) const ;
    void setTempFactor(float tf) { this->tempFactor = tf; };
};

std::vector<float3d> extractCoordinates(const std::vector<AtomCoordinate>& atoms);

static inline void extractCoordinates(
    float3d* output,
    const AtomCoordinate& atom1,
    const AtomCoordinate& atom2,
    const AtomCoordinate& atom3
) {
    output[0] = atom1.coordinate;
    output[1] = atom2.coordinate;
    output[2] = atom3.coordinate;
}

std::vector<AtomCoordinate> extractChain(
    std::vector<AtomCoordinate>& atoms, std::string chain
);

std::vector<AtomCoordinate> filterBackbone(const tcb::span<AtomCoordinate>& atoms);

void printAtomCoordinateVector(std::vector<AtomCoordinate>& atoms, int option = 0);

std::vector<AtomCoordinate> weightedAverage(
    const std::vector<AtomCoordinate>& origAtoms, const std::vector<AtomCoordinate>& revAtoms
);

// Length of one formatted PDB ATOM record (fixed columns + trailing '\n').
static constexpr int PDB_ATOM_LINE_LEN = 81;

// When `precomputedLines` is non-null it must point at an array of
// PDB_ATOM_LINE_LEN-byte ATOM lines (one per atom, in `atoms` order); the ATOM
// record for each atom is copied from there instead of being formatted here, and
// `*lineCursor` is advanced. Used by the GPU PDB writer. Default (nullptr)
// formats on the host exactly as before.
void writeAtomCoordinatesToPDB(
    const std::vector<AtomCoordinate>& atoms, const std::string& title, std::string& output,
    bool appendOutput = false, bool emitFinalTer = true,
    const char* precomputedLines = nullptr, size_t* lineCursor = nullptr
);
int writeAtomCoordinatesToPDBFile(
    const std::vector<AtomCoordinate>& atoms, const std::string& title, const std::string& pdb_path
);

#ifdef FOLDCOMP_WITH_MMCIF_OUTPUT
bool writeAtomCoordinatesToMMCIF(
    std::vector<AtomCoordinate>& atoms, const std::string& title, std::string& output
);
#endif

std::vector<std::vector<AtomCoordinate>> splitAtomByResidue(
    const tcb::span<AtomCoordinate>& atomCoordinates
);

std::vector<std::pair<size_t, size_t>> splitResidueRanges(
    const tcb::span<AtomCoordinate>& atomCoordinates
);

std::vector<std::string> getResidueNameVector(
    const tcb::span<AtomCoordinate>& atomCoordinates
);

AtomCoordinate findFirstAtom(const std::vector<AtomCoordinate>& atoms, std::string atom_name);
AtomCoordinate findFirstAtom(const tcb::span<const AtomCoordinate>& atoms, std::string atom_name);
void setAtomIndexSequentially(std::vector<AtomCoordinate>& atoms, int start);
void removeAlternativePosition(std::vector<AtomCoordinate>& atoms);
bool startsNewResidue(const AtomCoordinate& current, const AtomCoordinate& previous);

std::vector<AtomCoordinate> getAtomsWithResidueIndex(
    const tcb::span<AtomCoordinate>& atoms, int residue_index,
    std::vector<std::string> atomNames = {"N", "CA", "C"}
);

std::vector<std::vector<AtomCoordinate>> getAtomsWithResidueIndex(
    const tcb::span<AtomCoordinate>& atoms, std::vector<int> residue_index,
    std::vector<std::string> atomNames = {"N", "CA", "C"}
);
float RMSD(std::vector<AtomCoordinate>& atoms1, std::vector<AtomCoordinate>& atoms2);

template <int32_t T, int32_t P>
void ftoa(float n, char* s);

std::vector<std::pair<size_t, size_t>> identifyChains(const std::vector<AtomCoordinate>& atoms);
std::vector<std::pair<size_t, size_t>> identifyDiscontinousResInd(const std::vector<AtomCoordinate>& atoms, size_t chain_start, size_t chain_end);
std::vector<std::pair<size_t, size_t>> identifyCompleteBackboneRegions(const tcb::span<AtomCoordinate>& atoms);

struct BackboneRegion {
    size_t start;
    size_t end;
    bool encodable;
};

std::vector<BackboneRegion> identifyBackboneRegions(const tcb::span<AtomCoordinate>& atoms);

bool serializeAtomCoordinates(
    const std::vector<AtomCoordinate>& atoms,
    std::string& output
);

bool serializeAtomCoordinates(
    const tcb::span<const AtomCoordinate>& atoms,
    std::string& output
);

bool deserializeAtomCoordinates(
    const char* data,
    size_t size,
    std::vector<AtomCoordinate>& atoms
);
