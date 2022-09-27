/**
 * File: atom_coordinate.h
 * Project: foldcomp
 * Created: 2021-01-18 12:43:08
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     The data type to handle atom coordinate comes here.
 * ---
 * Last Modified: 2022-09-27 11:50:15
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */

#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

class AtomCoordinate {
private:
    void check3dCoordinate();
public:
    AtomCoordinate();
    AtomCoordinate(
        std::string a, std::string r, std::string c,
        int ai, int ri, float x, float y, float z
    );
    AtomCoordinate(
        std::string a, std::string r, std::string c,
        int ai, int ri, std::vector<float> coord
    );
    ~AtomCoordinate();
    // data
    std::string atom;
    std::string residue;
    std::string chain;
    int atom_index;
    int residue_index;
    std::vector<float> coordinate;
    float occupancy;
    float tempFactor;

    // operators
    bool operator==(const AtomCoordinate& other) const;
    bool operator!=(const AtomCoordinate& other) const;
    //BackboneChain toCompressedResidue();

    //methods
    bool isBackbone();
    void print(int option = 0);
    void setTempFactor(float tf) { this->tempFactor = tf; };
};

std::vector< std::vector<float> > extractCoordinates(std::vector<AtomCoordinate>& atoms);
std::vector<AtomCoordinate> extractChain(
    std::vector<AtomCoordinate>& atoms, std::string chain
);

std::vector<AtomCoordinate> filterBackbone(std::vector<AtomCoordinate>& atoms);

void printAtomCoordinateVector(std::vector<AtomCoordinate>& atoms, int option = 0);

std::vector<AtomCoordinate> weightedAverage(
    std::vector<AtomCoordinate>& origAtoms, std::vector<AtomCoordinate>& revAtoms
);

void writeAtomCoordinatesToPDB(
    std::vector<AtomCoordinate>& atoms, std::string title, std::ostream& pdb_path
);
int writeAtomCoordinatesToPDBFile(
    std::vector<AtomCoordinate>& atoms, std::string title, std::string pdb_path
);

std::vector< std::vector<AtomCoordinate> > splitAtomByResidue(
    const std::vector<AtomCoordinate>& atomCoordinates
);

std::vector<std::string> getResidueNameVector(
    const std::vector<AtomCoordinate>& atomCoordinates
);

AtomCoordinate findFirstAtom(std::vector<AtomCoordinate>& atoms, std::string atom_name);
void setAtomIndexSequentially(std::vector<AtomCoordinate>& atoms, int start);
void removeAlternativePosition(std::vector<AtomCoordinate>& atoms);

std::vector<AtomCoordinate> getAtomsWithResidueIndex(
    std::vector<AtomCoordinate>& atoms, int residue_index,
    std::vector<std::string> atomNames = {"N", "CA", "C"}
);

std::vector< std::vector<AtomCoordinate> > getAtomsWithResidueIndex(
    std::vector<AtomCoordinate>& atoms, std::vector<int> residue_index,
    std::vector<std::string> atomNames = {"N", "CA", "C"}
);

float distance(AtomCoordinate& a, AtomCoordinate& b);
std::vector<float> distance(std::vector<AtomCoordinate>& a, std::vector<AtomCoordinate>& b);
float RMSD(std::vector<AtomCoordinate>& atoms1, std::vector<AtomCoordinate>& atoms2);