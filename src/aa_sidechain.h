/**
 * File: aa_sidechain.h
 * Project: foldcomp
 * Created: 2021-07-07 13:38:27
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     Functions for handling side-chain atoms of amino acids
 * ---
 * Last Modified: 2022-07-20 02:00:18
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */

#pragma once
#include <string>
#include <vector>
#include <map>
#include "amino_acid.h"
#include "atom_coordinate.h"
#include "nerf.h"
#include "torsion.xyz.h"
#include "utility.h"

/**
 * @brief Reconstruct side chain
 * @param backbone
 * @return std::vector<AtomCoordinate>
 */
std::vector<AtomCoordinate> reconstructSideChain(std::vector<AtomCoordinate> backbone, AtomCoordinate c0);
std::vector<AtomCoordinate> reconstructSideChain(
    std::vector< std::vector<float> > backboneCoordinates,
    std::string residueName
);

// A Temporary function to check if BA, BL, TA calculation works
// 2021-09-30 20:57:36 Currently, this doesn't work either
std::vector<AtomCoordinate> reconstructSideChainFromCalculatedInfo(
    std::vector<AtomCoordinate> backbone
);

std::map<std::string, float> calculateBondLengths(
    std::vector<AtomCoordinate> originalAtoms, AminoAcid AA
);
std::map<std::string, float> calculateBondAngles(
    std::vector<AtomCoordinate> originalAtoms, AminoAcid AA
);
std::map<std::string, float> calculateTorsionAngles(
    std::vector<AtomCoordinate> originalAtoms, AminoAcid AA
);

std::vector<float> calculateTorsionAnglesInResidue(
    std::vector<AtomCoordinate>& originalAtoms, AminoAcid& AA
);

std::vector< std::vector<float> > calculateSideChainTorsionAnglesPerResidue(
    std::vector<AtomCoordinate>& originalAtoms, std::map<std::string, AminoAcid>& AAmap
);

std::map<std::string, AtomCoordinate> convertAtomCoordinateVectorToMap(
    std::vector<AtomCoordinate> input
);
void compareMap(
    std::map<std::string, float> m1, std::map<std::string, float> m2,
    std::string m1Name, std::string m2Name
);

// Trying to use pointers for the parameter
// "Pass by reference"

/**
 * @brief Calculate all the side-chain information
 *
 * @param peptideCoordinates
 * @param aminoAcid
 * @return int
 */
int calculateSideChainInfo(
    std::vector<AtomCoordinate>& peptideCoordinates,
    std::map<std::string, AminoAcid>& AAS
);

int reconstructSideChainFromPeptide(
    std::vector<AtomCoordinate>& peptideCoordinates,
    std::map<std::string, AminoAcid>& AAS,
    std::vector<AtomCoordinate>& output
);

std::vector<AtomCoordinate> residueReconstruction(
    std::vector<AtomCoordinate>& currResidue,
    AminoAcid& currAA,
    std::map<std::string, float>& currResTorsionAngles
);

std::vector<AtomCoordinate> residueReconstruction(
    std::vector<AtomCoordinate>& currResidue,
    std::string residueName,
    std::vector<float>& currResTorsionAngles
);

int saveTorsionAngles(
    std::vector<AtomCoordinate>& peptideCoordinates,
    std::map<std::string, AminoAcid>& AAS
);

void checkEmptyAtomsInResidue(std::vector<AtomCoordinate>& Residue, AminoAcid& AA);

std::map<std::string, std::vector< std::vector<float> > > groupSideChainTorsionByResidue(
    std::vector< std::vector<float> >& sideChainTorsionAngles,
    std::vector<std::string>& residueNames,
    std::map<std::string, AminoAcid>& AAS
);
std::vector<float> getSpecificTorsionAngle(
    std::map<std::string, std::vector< std::vector<float> > > sideChainTorsionMap,
    std::string& residueName, int index
);