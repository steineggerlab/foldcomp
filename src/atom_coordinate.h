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
#include <string>
#include <vector>

class AtomCoordinate {
public:
    AtomCoordinate() = default;
    AtomCoordinate(
        const std::string& a, const std::string& r, const std::string& c,
        int ai, int ri, float x, float y, float z,
        float occupancy = 0.0f, float tempFactor = 0.0f,
        int model = 1, char insertionCode = ' ', char altloc = ' '
    );
    AtomCoordinate(
        const std::string& a, const std::string& r, const std::string& c,
        int ai, int ri, float3d coord,
        float occupancy = 0.0f, float tempFactor = 0.0f,
        int model = 1, char insertionCode = ' ', char altloc = ' '
    );
    // data
    std::string atom;
    std::string residue;
    std::string chain;
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

void writeAtomCoordinatesToPDB(
    const std::vector<AtomCoordinate>& atoms, const std::string& title, std::string& output,
    bool appendOutput = false, bool emitFinalTer = true
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
