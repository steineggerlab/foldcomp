/**
 * File: foldcomp.cpp
 * Project: foldcomp
 * Created: 2021-02-04 13:31:52
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     This file contains main data structures for torsion angle compression and
 *     functions for handling them.
 * ---
 * Last Modified: 2025-04-21 15:37:18
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
#include "foldcomp.h"

#include "amino_acid.h"
#include "sidechain.h"
#include "torsion_angle.h"
#include "utility.h"

#include <cmath>

#include <algorithm>
#include <bitset>
#include <cstdio>
#include <utility>

namespace {

template <typename T>
void appendBytes(std::string& output, const T& value) {
    size_t offset = output.size();
    output.resize(offset + sizeof(T));
    memcpy(output.data() + offset, &value, sizeof(T));
}

template <typename T>
void appendBytes(std::string& output, const std::vector<T>& values) {
    if (values.empty()) {
        return;
    }
    size_t bytes = values.size() * sizeof(T);
    size_t offset = output.size();
    output.resize(offset + bytes);
    memcpy(output.data() + offset, values.data(), bytes);
}

template <typename T>
void appendByteVector(std::string& output, const std::vector<T>& values) {
    if (values.empty()) {
        return;
    }
    size_t offset = output.size();
    output.resize(offset + values.size());
    for (size_t i = 0; i < values.size(); i++) {
        output[offset + i] = static_cast<char>(values[i]);
    }
}

}

// Changed at 2022-04-14 16:02:30
/**
 * @brief Convert a BackboneChain to a byte array (8 bytes)
 *
 * @param res a BackboneChain
 * @return char*
 */
int convertBackboneChainToBytes(BackboneChain& res, char* output) {
    int flag = 0;
    // 00-04: residue, 05-07: OMEGA first 3 bits
    output[0] = ((res.residue << 3) | ((res.omega & 0x07FF) >> 8));
    // 08-15: OMEGA last 8 bits
    output[1] = res.omega & 0x00FF;
    // 16-23: PSI first 8 bits
    output[2] = ((res.psi & 0x0FFF) >> 4);
    // 24-27: PSI last 4 bits, 28-31: PHI first 4 bits
    output[3] = ((res.psi & 0x000F) << 4) | ((res.phi & 0x0FFF) >> 8);
    // 32-39: PHI last 8 bits
    output[4] = res.phi & 0x00FF;
    // 40-47: CA_C_N
    output[5] = res.ca_c_n_angle;
    // 48-55: C_N_CA
    output[6] = res.c_n_ca_angle;
    // 56-63: N_CA_C
    output[7] = res.n_ca_c_angle;
    return flag;
}

/**
 * @brief Read a byte array and convert it to a BackboneChain.
 * This function is used for reading compressed residue from a file.
 * @param bytes a byte array (8 bytes) which encodes a BackboneChain
 * @return BackboneChain
 */
BackboneChain convertBytesToBackboneChain(char* bytes) {
    BackboneChain res;
    // 00-04: residue
    res.residue = ((bytes[0] & 0xF8) >> 3);
    // 05-07: OMEGA first 3 bits, 08-15: OMEGA last 8 bits
    res.omega = (unsigned int)(((bytes[0] & 0x0007) << 8) | (bytes[1] & 0x00FF));
    // 16-23: PSI first 8 bits, 24-27: PSI last 4 bits
    res.psi = (unsigned int)(((bytes[2] & 0x00FF) << 4) | (bytes[3] & 0x00FF) >> 4);
    // 28-31: PHI first 4 bits, 32-39: PHI last 8 bits
    res.phi = (unsigned int)(((bytes[3] & 0x000F) << 8) | (bytes[4] & 0x00FF));
    // 40-47: CA_C_N
    res.ca_c_n_angle = bytes[5];
    // 48-55: C_N_CA
    res.c_n_ca_angle = bytes[6];
    // 56-63: N_CA_C
    res.n_ca_c_angle = bytes[7];
    return res;
}

// NOTE:
BackboneChain newBackboneChain(
    char residue, unsigned int phi, unsigned int psi, unsigned int omega,
    unsigned int n_ca_c_angle, unsigned int ca_c_n_angle, unsigned int c_n_ca_angle
) {
    // unsigned int r = convertOneLetterCodeToInt(residue);
    BackboneChain res;
    res.residue = convertOneLetterCodeToInt(residue);
    res.ca_c_n_angle = ca_c_n_angle;
    res.c_n_ca_angle = c_n_ca_angle;
    res.n_ca_c_angle = n_ca_c_angle;
    res.psi = psi;
    res.omega = omega;
    res.phi = phi;
    return res;
}

BackboneChain newBackboneChain(
    unsigned int bResidue, unsigned int phi, unsigned int psi, unsigned int omega,
    unsigned int n_ca_c_angle, unsigned int ca_c_n_angle, unsigned int c_n_ca_angle
) {
    BackboneChain res;
    res.residue = bResidue;
    res.ca_c_n_angle = ca_c_n_angle;
    res.c_n_ca_angle = c_n_ca_angle;
    res.n_ca_c_angle = n_ca_c_angle;
    res.psi = psi;
    res.omega = omega;
    res.phi = phi;
    return res;
}


// WARNING:
/**
 * @brief Convert a BackboneChain to DecompressedBackboneChain
 *
 * @description As the angles are short-encoded in the compressed format,
 * this function converts the short-encoded angles to the float.
 * @param bb
 * @param header
 * @return DecompressedBackboneChain
 */
DecompressedBackboneChain decompressBackboneChain(
    const BackboneChain& bb, const CompressedFileHeader& header
) {
    DecompressedBackboneChain output;
    output.residue = convertIntToOneLetterCode(bb.residue);
    output.phi = _continuize(bb.phi, header.mins[0], header.cont_fs[0]);
    output.psi = _continuize(bb.psi, header.mins[1], header.cont_fs[1]);
    output.omega = _continuize(bb.omega, header.mins[2], header.cont_fs[2]);
    output.n_ca_c_angle = _continuize(bb.n_ca_c_angle, header.mins[3], header.cont_fs[3]);
    output.ca_c_n_angle = _continuize(bb.ca_c_n_angle, header.mins[4], header.cont_fs[4]);
    output.c_n_ca_angle = _continuize(bb.c_n_ca_angle, header.mins[5], header.cont_fs[5]);
    return output;
}


/**
 * @brief Convert a vectorof BackboneChain to DecompressedBackboneChain vector
 *
 * @param bbv
 * @param header
 * @return std::vector<DecompressedBackboneChain>
 */
std::vector<DecompressedBackboneChain> decompressBackboneChain(
    const std::vector<BackboneChain>& bbv, const CompressedFileHeader& header
) {
    std::vector<DecompressedBackboneChain> output;
    output.reserve(bbv.size());
    for (const auto& bb : bbv) {
        output.push_back(decompressBackboneChain(bb, header));
    }
    return output;
}

float _continuize(unsigned int input, float min, float cont_f) {
    float output = min + ((float)input * cont_f);
    return output;
}

std::string getChainName(const CompressedFileHeader& header) {
    std::string chain;
    if (header.chain != '\0') {
        chain.push_back(header.chain);
    }
    if (header.chain2 != '\0') {
        chain.push_back(header.chain2);
    }
    if (header.chain3 != '\0') {
        chain.push_back(header.chain3);
    }
    return chain;
}

void setChainName(CompressedFileHeader& header, const std::string& chain) {
    header.chain = chain.empty() ? '\0' : chain[0];
    header.chain2 = chain.size() > 1 ? chain[1] : '\0';
    header.chain3 = chain.size() > 2 ? chain[2] : '\0';
}

/**
 * @brief Reconstruct backbone atoms from compressed backbone info
 *
 * @param prev_atoms std::vector<AtomCoordinates> of 3 previous atoms
 * @param backbone std::vector<BackboneChain>
 * @return std::vector<AtomCoordinate>
 */
std::vector<AtomCoordinate> reconstructBackboneAtoms(
    const std::vector<AtomCoordinate>& prevAtoms,
    const std::vector<BackboneChain>& backbone,
    CompressedFileHeader& header
) {
    Nerf nerf;
    // Save first three atoms
    std::vector<AtomCoordinate> reconstructedAtoms = {
        prevAtoms[0], prevAtoms[1], prevAtoms[2]
    };
    int total = backbone.size();

    std::vector<DecompressedBackboneChain> deBackbone = decompressBackboneChain(backbone, header);
    int currAtomIndex = prevAtoms[2].atom_index + 1;
    int currResidueIndex = prevAtoms[2].residue_index + 1;

    // Iterate through backbone
    // Should put N, CA, C in this loop
    for (int i = 0; i < (total - 1); i++) {
        const AtomCoordinate& prevAtom1 = reconstructedAtoms[i*3];
        const AtomCoordinate& prevAtom2 = reconstructedAtoms[i*3 + 1];
        const AtomCoordinate& prevAtom3 = reconstructedAtoms[i*3 + 2];
        float3d prevCoords[3];
        extractCoordinates(prevCoords, prevAtom1, prevAtom2, prevAtom3);
        // Convert char (deBackbone[i].residue) to string (currResidue)
        std::string currResidue = getThreeLetterCode(deBackbone[i + 1].residue);
        std::string currChain = prevAtom1.chain;

        // Place N
        float3d currNCoord = nerf.place_atom(
            prevCoords, C_TO_N_DIST, deBackbone[i].ca_c_n_angle, deBackbone[i].psi
        );
        // Place CA
        prevCoords[0] = prevCoords[1];
        prevCoords[1] = prevCoords[2];
        prevCoords[2] = currNCoord;
        float3d currCACoord;
        if (deBackbone[i].residue != 'P') {
            currCACoord = nerf.place_atom(
                prevCoords, N_TO_CA_DIST, deBackbone[i].c_n_ca_angle, deBackbone[i].omega
            );
        } else {
            currCACoord = nerf.place_atom(
                prevCoords, PRO_N_TO_CA_DIST, deBackbone[i].c_n_ca_angle, deBackbone[i].omega
            );
        }

        // Place C
        currAtomIndex++;
        prevCoords[0] = prevCoords[1];
        prevCoords[1] = prevCoords[2];
        prevCoords[2] = currCACoord;
        float3d currCCoord = nerf.place_atom(
            prevCoords, CA_TO_C_DIST, deBackbone[i].n_ca_c_angle, deBackbone[i].phi
        );

        reconstructedAtoms.emplace_back(
            "N", currResidue, currChain,
            currAtomIndex, currResidueIndex,
            currNCoord.x, currNCoord.y, currNCoord.z
        );
        currAtomIndex++;
        reconstructedAtoms.emplace_back(
            "CA", currResidue, currChain,
            currAtomIndex, currResidueIndex,
            currCACoord.x, currCACoord.y, currCACoord.z
        );
        currAtomIndex++;
        reconstructedAtoms.emplace_back(
            "C", currResidue, currChain,
            currAtomIndex, currResidueIndex,
            currCCoord.x, currCCoord.y, currCCoord.z
        );
        // Increment Residue Index
        currResidueIndex++;
        currAtomIndex++;

    }
    return reconstructedAtoms;
}

int reconstructBackboneReverse(
    std::vector<AtomCoordinate>& atom, std::vector< std::vector<float> >& lastCoords,
    std::vector<float>& torsion_angles, Nerf& nerf
) {
    std::vector<AtomCoordinate> forwardAtom = atom;
    // Last atoms
    atom[atom.size() - 3].coordinate.x = lastCoords[0][0];
    atom[atom.size() - 3].coordinate.y = lastCoords[0][1];
    atom[atom.size() - 3].coordinate.z = lastCoords[0][2];
    atom[atom.size() - 2].coordinate.x = lastCoords[1][0];
    atom[atom.size() - 2].coordinate.y = lastCoords[1][1];
    atom[atom.size() - 2].coordinate.z = lastCoords[1][2];
    atom[atom.size() - 1].coordinate.x = lastCoords[2][0];
    atom[atom.size() - 1].coordinate.y = lastCoords[2][1];
    atom[atom.size() - 1].coordinate.z = lastCoords[2][2];

    std::vector<float> bond_angles = nerf.getBondAngles(forwardAtom);

    std::vector<AtomCoordinate> atomBackward = nerf.reconstructWithReversed(
        atom, torsion_angles, bond_angles
    );

    int total = atom.size();
    for (int i = 0; i < total; i++) {
        atom[i].coordinate.x =
            ((forwardAtom[i].coordinate.x * static_cast<float>(total - i)) +
             (atomBackward[i].coordinate.x * static_cast<float>(i))) / static_cast<float>(total);
        atom[i].coordinate.y =
            ((forwardAtom[i].coordinate.y * static_cast<float>(total - i)) +
             (atomBackward[i].coordinate.y * static_cast<float>(i))) / static_cast<float>(total);
        atom[i].coordinate.z =
            ((forwardAtom[i].coordinate.z * static_cast<float>(total - i)) +
             (atomBackward[i].coordinate.z * static_cast<float>(i))) / static_cast<float>(total);
    }
    return 0;
}


int discretizeSideChainTorsionAngles(
    std::vector< std::vector<float> >& torsionPerResidue,
    std::vector<std::string>& residueNames,
    const std::map<std::string, AminoAcid>& AAS,
    SideChainDiscretizers& scDiscretizers,
    std::map<std::string, std::vector<Discretizer> >& scDiscretizersMap,
    std::vector<unsigned int>& output
) {
    // Declare
    int success = 0;
    std::string currResidue;
    int currResidueTorsionNum;
    float min, cont_f;
    std::vector<float> currTorsion;
    std::vector<unsigned int> currTorsionDiscretized;
    float* min_arr;
    float* cont_f_arr;

    std::map<std::string, std::vector< std::vector<float> > > sideChainTorsionMap;
    sideChainTorsionMap = groupSideChainTorsionByResidue(torsionPerResidue, residueNames, AAS);

    // Fill in Discretizer map
    for (const auto& sc : sideChainTorsionMap) {
        currResidue = sc.first;
        currResidueTorsionNum = getSideChainTorsionNum(currResidue);
        min_arr = getMinPointerFromSideChainDiscretizers(currResidue, scDiscretizers);
        cont_f_arr = getContFFromSideChainDiscretizers(currResidue, scDiscretizers);
        for (int i = 0; i < currResidueTorsionNum; i++) {
            // Get current torsion angle vector
            currTorsion = getSpecificTorsionAngle(sideChainTorsionMap, currResidue, i);
            // Discretize the torsion angles
            // Discretizer currDiscretizer = Discretizer(currTorsion, pow(2, NUM_BITS_SIDECHAIN) - 1);
            // TODO: TESTING FIXEDANGLEDISCRETIZER
            FixedAngleDiscretizer currDiscretizer = FixedAngleDiscretizer(pow(2, NUM_BITS_TEMP) - 1);

            // Save min and cont_f to SideChainDiscretizers
            min = currDiscretizer.min;
            cont_f = currDiscretizer.cont_f;
            min_arr[i] = min;
            cont_f_arr[i] = cont_f;
            // Save Discretizer to map
            scDiscretizersMap[currResidue][i] = currDiscretizer;
        }
    }

    Discretizer torsionDisc;
    unsigned int torsionDiscretized;
    // Discretize the torsion angles and make the result into a flattend vector
    for (size_t i = 0; i < torsionPerResidue.size(); i++) {
        currResidue = residueNames[i];
        currResidueTorsionNum = getSideChainTorsionNum(currResidue);
        for (int j = 0; j < currResidueTorsionNum; j++) {
            torsionDisc = scDiscretizersMap[currResidue][j];
            torsionDiscretized = torsionDisc.discretize(torsionPerResidue[i][j]);
            // Append to flattened vector
            output.push_back(torsionDiscretized);
        }
    }

    return success;
}

int continuizeSideChainTorsionAngles(
    std::vector<unsigned int>& torsionDiscretized,
    std::vector<std::string>& residueNames,
    SideChainDiscretizers& scDiscretizers,
    std::map<std::string, std::vector<Discretizer> >& scDiscretizersMap,
    std::vector< std::vector<float> >& output
) {

    scDiscretizersMap = initializeSideChainDiscMap();
    int success = fillSideChainDiscretizerMap(scDiscretizers, scDiscretizersMap);
    std::vector< std::vector<float> > torsionPerResidue;
    std::vector<float> currTorsionVector;
    FixedAngleDiscretizer currDiscretizer(pow(2, NUM_BITS_TEMP) - 1);
    // Iterate
    int currIndex = 0;
    for (size_t i = 0; i < residueNames.size(); i++) {
        std::string currResidue = residueNames[i];
        int currResidueTorsionNum = getSideChainTorsionNum(currResidue);
        currTorsionVector.clear();
        currTorsionVector.resize(currResidueTorsionNum);
        for (int j = 0; j < currResidueTorsionNum; j++) {
            // Get current torsion angle vector
            //currTorsion = scDiscretizersMap[currResidue][j].continuize(torsionDiscretized[currIndex]);
            float currTorsion = currDiscretizer.continuize(torsionDiscretized[currIndex]);
            currTorsionVector[j] = currTorsion;
            currIndex++;
        }
        torsionPerResidue.push_back(currTorsionVector);
    }
    output = torsionPerResidue;
    return success;
}

int fillSideChainDiscretizerMap(
    SideChainDiscretizers& scDiscretizers,
    std::map<std::string, std::vector<Discretizer> >& scDiscretizersMap
) {
    // Declare
    int success = 0;
    std::string currResidue;
    int currResidueTorsionNum;
    float min, cont_f;
    // Iterate through map
    for (const auto& sc : scDiscretizersMap) {
        currResidue = sc.first;
        currResidueTorsionNum = getSideChainTorsionNum(currResidue);
        for (int i = 0; i < currResidueTorsionNum; i++) {
            // Get current torsion angle vector
            min = getMinPointerFromSideChainDiscretizers(currResidue, scDiscretizers)[i];
            cont_f = getContFFromSideChainDiscretizers(currResidue, scDiscretizers)[i];
            // Save min and cont_f to SideChainDiscretizers
            scDiscretizersMap[currResidue][i].min = min;
            scDiscretizersMap[currResidue][i].cont_f = cont_f;
        }
    }
    return success;
}

int Foldcomp::_discretizeSideChainTorsionAngles(
    std::vector< std::vector<float> >& input,
    std::vector<unsigned int>& output
) {
    // Declare
    int success = 0;
    success = discretizeSideChainTorsionAngles(
        input, this->residueThreeLetter, this->AAS,
        this->sideChainDisc, this->sideChainDiscMap, output
    );
    this->header.nSideChainTorsion = output.size();
    this->nSideChainTorsion = output.size();
    return success;
}

int Foldcomp::_continuizeSideChainTorsionAngles(
    std::vector<unsigned int>& input, std::vector< std::vector<float> >& output
) {
    // Declare
    int success = 0;
    success = continuizeSideChainTorsionAngles(
        input, this->residueThreeLetter, this->sideChainDisc, this->sideChainDiscMap, output
    );
    return success;
}

/**
 * @brief [TEMP] Print the bit representation of the variables in the
 * compressed residue
 *
 * @param res
 */
void printCompressedResidue(BackboneChain& res) {
    // Print the variables in the header
    std::cout << "SIZE: " << sizeof(res) << std::endl;
    std::cout << "residue: " << res.residue << std::endl;
    std::cout << "phi: " << res.phi << std::endl;
    std::cout << "psi: " << res.psi << std::endl;
    std::cout << "omega: " << res.omega << std::endl;
    std::cout << "n_ca_c_angle: " << res.n_ca_c_angle << std::endl;
    std::cout << "ca_c_n_angle: " << res.ca_c_n_angle << std::endl;
    std::cout << "c_n_ca_angle: " << res.c_n_ca_angle << std::endl;
    std::cout << "CONVERTED BYTE ARRAY: ";
    std::bitset<8> bits;
    char* byteArray = new char[8];
    // int flag = convertBackboneChainToBytes(res, byteArray);
    for (int i = 0; i < 8; i++) {
        bits = byteArray[i];
        std::cout << bits << " ";
    }
    delete[] byteArray;
    std::cout << std::endl;
}

int Foldcomp::preprocess(const tcb::span<AtomCoordinate>& atoms) {
    int success = 0;
    if (atoms.empty()) {
        return -1;
    }
    this->isPreprocessed = false;
    this->isCompressed = false;
    this->compressedBackBone.clear();
    this->compressedSideChain.clear();
    this->residues.clear();
    this->backbone.clear();
    this->residueThreeLetter.clear();
    this->backboneTorsionAngles.clear();
    this->backboneBondAngles.clear();
    this->psi.clear();
    this->omega.clear();
    this->phi.clear();
    this->n_ca_c_angle.clear();
    this->ca_c_n_angle.clear();
    this->c_n_ca_angle.clear();
    this->psiDiscretized.clear();
    this->omegaDiscretized.clear();
    this->phiDiscretized.clear();
    this->n_ca_c_angleDiscretized.clear();
    this->ca_c_n_angleDiscretized.clear();
    this->c_n_ca_angleDiscretized.clear();
    this->sideChainAnglesPerResidue.clear();
    this->sideChainAnglesDiscretized.clear();
    this->tempFactors.clear();
    this->tempFactorsDiscretized.clear();
    this->anchorAtoms.clear();
    this->anchorCoordinates.clear();
    this->anchorIndices.clear();

    // Discretize
    // Keep only residues with complete backbone (N/CA/C).
    std::vector<AtomCoordinate> encodableAtoms;
    std::vector<std::pair<size_t, size_t>> residueRanges = splitResidueRanges(atoms);
    encodableAtoms.reserve(atoms.size());
    this->backbone.reserve(residueRanges.size() * 3);
    this->residueThreeLetter.reserve(residueRanges.size());
    this->tempFactors.reserve(residueRanges.size());
    for (const auto& residue : residueRanges) {
        const AtomCoordinate* n = nullptr;
        const AtomCoordinate* ca = nullptr;
        const AtomCoordinate* c = nullptr;
        for (size_t i = residue.first; i < residue.second; i++) {
            const auto& atom = atoms[i];
            if (atom.atom == "N") {
                if (n == nullptr) {
                    n = &atom;
                }
            } else if (atom.atom == "CA") {
                if (ca == nullptr) {
                    ca = &atom;
                }
            } else if (atom.atom == "C") {
                if (c == nullptr) {
                    c = &atom;
                }
            }
        }
        if (n != nullptr && ca != nullptr && c != nullptr) {
            encodableAtoms.insert(encodableAtoms.end(), atoms.begin() + residue.first, atoms.begin() + residue.second);
            this->backbone.emplace_back(*n);
            this->backbone.emplace_back(*ca);
            this->backbone.emplace_back(*c);
            this->residueThreeLetter.push_back(atoms[residue.first].residue);
            this->tempFactors.push_back(ca->tempFactor);
        }
    }
    if (encodableAtoms.empty()) {
        return -1;
    }
    if (this->backbone.size() < 6 || (this->backbone.size() % 3) != 0) {
        return -1;
    }
    for (size_t i = 0; i < this->backbone.size(); i += 3) {
        if (this->backbone[i].atom != "N" ||
            this->backbone[i + 1].atom != "CA" ||
            this->backbone[i + 2].atom != "C") {
            return -1;
        }
    }

    // this->title = this->strTitle.c_str();
    this->lenTitle = this->strTitle.size();

    this->nResidue = this->backbone.size() / 3;
    this->nBackbone = this->backbone.size();
    this->nAtom = encodableAtoms.size();
    this->idxResidue = this->backbone[0].residue_index;
    this->idxAtom = this->backbone[0].atom_index;
    this->chain = this->backbone[0].chain.c_str()[0];
    this->firstResidue = getOneLetterCode(this->backbone[0].residue);
    this->lastResidue = getOneLetterCode(this->backbone[this->backbone.size() - 1].residue);

    // Anchor atoms
    this->_setAnchor();

    if (encodableAtoms[encodableAtoms.size() - 1].atom == "OXT") {
        this->hasOXT = 1;
        this->OXT = encodableAtoms[encodableAtoms.size() - 1];
        this->OXT_coords = encodableAtoms[encodableAtoms.size() - 1].coordinate;
    } else {
        this->hasOXT = 0;
        this->OXT = AtomCoordinate();
        this->OXT_coords = {0.0, 0.0, 0.0};
    }

    std::vector<float> backboneTorsion = getTorsionFromXYZ(this->backbone, 1);
    this->backboneTorsionAngles = backboneTorsion;
    // Split backbone into phi, psi, omega
    // Calculate phi, psi, omega
    for (size_t i = 0; i < backboneTorsion.size(); i += 3) {
        this->psi.push_back(backboneTorsion[i]);
        this->omega.push_back(backboneTorsion[i + 1]);
        this->phi.push_back(backboneTorsion[i + 2]);
    }

    std::vector<float> backboneBondAngles = this->nerf.getBondAngles(backbone);
    this->backboneBondAngles = backboneBondAngles;
    // Split bond angles into three parts
    for (size_t i = 1; i < backboneBondAngles.size(); i++) {
        if (i % 3 == 0) {
            this->n_ca_c_angle.push_back(backboneBondAngles[i]);
        } else if (i % 3 == 1) {
            this->ca_c_n_angle.push_back(backboneBondAngles[i]);
        } else {
            this->c_n_ca_angle.push_back(backboneBondAngles[i]);
        }
    }

    // Discretize
    this->phiDisc = Discretizer(this->phi, pow(2, NUM_BITS_PHI_PSI) - 1);
    this->phiDiscretized = this->phiDisc.discretize(this->phi);
    this->omegaDisc = FixedAngleDiscretizer(pow(2, NUM_BITS_OMEGA) - 1);
    this->omegaDiscretized = this->omegaDisc.discretize(this->omega);
    this->psiDisc = Discretizer(this->psi, pow(2, NUM_BITS_PHI_PSI) - 1);
    this->psiDiscretized = this->psiDisc.discretize(this->psi);
    this->n_ca_c_angleDisc = Discretizer(this->n_ca_c_angle, pow(2, NUM_BITS_BOND) - 1);
    this->n_ca_c_angleDiscretized = this->n_ca_c_angleDisc.discretize(this->n_ca_c_angle);
    this->ca_c_n_angleDisc = Discretizer(this->ca_c_n_angle, pow(2, NUM_BITS_BOND) - 1);
    this->ca_c_n_angleDiscretized = this->ca_c_n_angleDisc.discretize(this->ca_c_n_angle);
    this->c_n_ca_angleDisc = Discretizer(this->c_n_ca_angle, pow(2, NUM_BITS_BOND) - 1);
    this->c_n_ca_angleDiscretized = this->c_n_ca_angleDisc.discretize(this->c_n_ca_angle);

    // Set Discretizer for side chain
    this->sideChainDiscMap = initializeSideChainDiscMap();
    // Calculate side chain info

    this->sideChainAnglesPerResidue = calculateSideChainTorsionAnglesPerResidue(encodableAtoms, this->AAS);

    // Discretize side chain
    //this->_discretizeSideChainTorsionAngles(this->sideChainAnglesPerResidue, this->sideChainAnglesDiscretized);
    FixedAngleDiscretizer sideChainDiscretizer = FixedAngleDiscretizer(pow(2, NUM_BITS_TEMP) - 1);
    for (size_t i = 0; i < this->sideChainAnglesPerResidue.size(); i++) {
        for (size_t j = 0; j < this->sideChainAnglesPerResidue[i].size(); j++) {
            unsigned int temp = sideChainDiscretizer.discretize(this->sideChainAnglesPerResidue[i][j]);
            this->sideChainAnglesDiscretized.push_back(temp);
        }
    }
    this->nSideChainTorsion = this->sideChainAnglesDiscretized.size();

    // Get tempFactors
    // 2022-08-31 16:28:30 - Changed to save one tempFactor per residue
    // Discretize
    this->tempFactorsDisc = Discretizer(this->tempFactors, pow(2, NUM_BITS_TEMP) - 1);
    this->tempFactorsDiscretized = this->tempFactorsDisc.discretize(this->tempFactors);

    // Get header
    this->header = this->get_header();

    // Mark as processed
    this->isPreprocessed = true;

    return success;
}


std::vector<BackboneChain> Foldcomp::compress(const tcb::span<AtomCoordinate>& atoms) {
    std::vector<BackboneChain> output;
    this->residues.clear();
    // TODO: convert the atom coordinate vector into a vector of compressed residue
    // CURRENT VERSION - 2022-01-10 15:34:21
    // IGNORE BREAKS
    if (!this->isPreprocessed) {
        if (this->preprocess(atoms) != 0) {
            return output;
        }
    }
    if (atoms.size() < 3 || this->backbone.size() < 3 || this->nResidue <= 0) {
        return output;
    }
    output.reserve(static_cast<size_t>(this->nResidue));
    this->residues.reserve(static_cast<size_t>(this->nResidue));
    this->prevAtoms = {this->backbone[0], this->backbone[1], this->backbone[2]};

    AtomCoordinate currN;
    char currResCode;

    this->lastAtoms = {this->backbone[this->backbone.size() - 1],
                       this->backbone[this->backbone.size() - 2],
                       this->backbone[this->backbone.size() - 3]};

    // Need to extract backbone atoms in a separate vector

    BackboneChain res;
    for (int i = 0; i < (this->nResidue - 1); i++) {
        currN = this->backbone[i * 3];
        currResCode = getOneLetterCode(currN.residue);
        this->residues.push_back(currResCode);
        res.residue = convertOneLetterCodeToInt(currResCode);
        res.psi = this->psiDiscretized[i];
        res.omega = this->omegaDiscretized[i];
        res.phi = this->phiDiscretized[i];
        res.n_ca_c_angle = this->n_ca_c_angleDiscretized[i];
        res.ca_c_n_angle = this->ca_c_n_angleDiscretized[i];
        res.c_n_ca_angle = this->c_n_ca_angleDiscretized[i];
        output.push_back(res);
    }
    currN = this->backbone[(this->nResidue - 1) * 3];
    currResCode = getOneLetterCode(currN.residue);
    this->residues.push_back(currResCode);
    res.residue = convertOneLetterCodeToInt(currResCode);
    res.psi = 0; res.omega = 0; res.phi = 0;
    res.n_ca_c_angle = 0; res.ca_c_n_angle = 0; res.c_n_ca_angle = 0;
    output.push_back(res);
    this->compressedBackBone = output;

    this->isCompressed = true;
    return output;
}

int _restoreResidueNames(
    std::vector<BackboneChain>& compressedBackbone,
    CompressedFileHeader& /* header */,
    std::vector<char>& residueOneLetter,
    std::vector<std::string>& residueThreeLetter
) {
    int success = 0;
    std::string threeLetterCode;
    char oneLetterCode;
    residueThreeLetter.clear();
    for (size_t i = 0; i < compressedBackbone.size(); i++) {
        oneLetterCode = convertIntToOneLetterCode(compressedBackbone[i].residue);
        threeLetterCode = convertIntToThreeLetterCode(compressedBackbone[i].residue);
        residueOneLetter.push_back(oneLetterCode);
        residueThreeLetter.push_back(threeLetterCode);
    }
    return success;
}

int Foldcomp::_restoreDiscretizer(int angleType) {
    int success = 0;
    std::vector<unsigned int> temp;
    temp.reserve(this->compressedBackBone.size());
    for (const auto& cb : this->compressedBackBone) {
        switch (angleType) {
        case 0: // Phi
            temp.push_back(cb.phi);
            break;
        case 1:
            temp.push_back(cb.psi);
            break;
        case 2:
            temp.push_back(cb.omega);
            break;
        case 3:
            temp.push_back(cb.n_ca_c_angle);
            break;
        case 4:
            temp.push_back(cb.ca_c_n_angle);
            break;
        case 5:
            temp.push_back(cb.c_n_ca_angle);
            break;
        default:
            break;
        }
    }
    // Set
    switch (angleType) {
    case 0: // Phi
        this->phiDiscretized = temp;
        this->phiDisc.min = this->header.mins[0];
        this->phiDisc.cont_f = this->header.cont_fs[0];
        break;
    case 1: // Psi
        this->psiDiscretized = temp;
        this->psiDisc.min = this->header.mins[1];
        this->psiDisc.cont_f = this->header.cont_fs[1];
        break;
    case 2: // Omega
        this->omegaDiscretized = temp;
        this->omegaDisc.min = this->header.mins[2];
        this->omegaDisc.cont_f = this->header.cont_fs[2];
        break;
    case 3: // N-CA-C
        this->n_ca_c_angleDiscretized = temp;
        this->n_ca_c_angleDisc.min = this->header.mins[3];
        this->n_ca_c_angleDisc.cont_f = this->header.cont_fs[3];
        break;
    case 4: // CA-C-N
        this->ca_c_n_angleDiscretized = temp;
        this->ca_c_n_angleDisc.min = this->header.mins[4];
        this->ca_c_n_angleDisc.cont_f = this->header.cont_fs[4];
        break;
    case 5: // C-N-CA
        this->c_n_ca_angleDiscretized = temp;
        this->c_n_ca_angleDisc.min = this->header.mins[5];
        this->c_n_ca_angleDisc.cont_f = this->header.cont_fs[5];
        break;
    default:
        break;
    }

    return success;
}

// TODO: Write a function to restore prev atoms (AtomCoordinate) from coordinates (float)
// 2022-02-17 23:29:22
/**
 * @brief Restore previous atoms from the header & coordinates
 *
 * @param coords
 * @return int
 */
int Foldcomp::_restoreAtomCoordinate(float* coords) {
    int success = 0;
    std::string firstResidue = getThreeLetterCode(this->header.firstResidue);
    // TODO: Fix the chain alphabet according to the real data
    // 2022-03-04 14:30:58
    // IMPORTANT: WARNING: TODO: Atom indexing is not correct right now
    // SIDE CHAIN ATOMS SHOULD GE CONSIDERED WHEN INDEXING
    // convert char to string
    std::string chain = getChainName(this->header);
    AtomCoordinate prevN = AtomCoordinate(
        "N", firstResidue, chain, this->header.idxAtom, this->header.idxResidue,
        coords[0], coords[1], coords[2]
    );
    AtomCoordinate prevCA = AtomCoordinate(
        "CA", firstResidue, chain, this->header.idxAtom + 1, this->header.idxResidue,
        coords[3], coords[4], coords[5]
    );
    AtomCoordinate prevC = AtomCoordinate(
        "C", firstResidue, chain, this->header.idxAtom + 2, this->header.idxResidue,
        coords[6], coords[7], coords[8]
    );

    // Check if this->prevAtoms is empty
    if (this->prevAtoms.size() == 0) {
        this->prevAtoms.push_back(prevN);
        this->prevAtoms.push_back(prevCA);
        this->prevAtoms.push_back(prevC);
    } else {
        // If not empty, it will be updated
        this->prevAtoms[0] = prevN;
        this->prevAtoms[1] = prevCA;
        this->prevAtoms[2] = prevC;
    }

    return success;
}

int Foldcomp::_getAnchorNum(int threshold) {
    int nAnchor = 0;
    nAnchor = this->nResidue / threshold;
    return nAnchor;
}

void Foldcomp::_setAnchor() {
    this->nInnerAnchor = this->_getAnchorNum(this->anchorThreshold);
    this->nAllAnchor = this->nInnerAnchor + 2; // Start and end
    // Set the anchor points - residue index
    this->anchorIndices.clear();
    int interval = this->nResidue / (this->nAllAnchor - 1);
    for (int i = 0; i < this->nAllAnchor - 1; i++) {
        this->anchorIndices.push_back(i * interval);
    }
    this->anchorIndices.push_back(this->nResidue - 1);
    //
    this->anchorAtoms.clear();
    this->anchorAtoms.reserve(this->anchorIndices.size());
    for (size_t i = 0; i < this->anchorIndices.size(); i++) {
        size_t idx = static_cast<size_t>(this->anchorIndices[i]) * 3;
        if (idx + 2 < this->backbone.size()) {
            this->anchorAtoms.push_back({
                this->backbone[idx],
                this->backbone[idx + 1],
                this->backbone[idx + 2]
            });
        } else if (this->backbone.size() >= 3) {
            this->anchorAtoms.push_back({
                this->backbone[this->backbone.size() - 3],
                this->backbone[this->backbone.size() - 2],
                this->backbone[this->backbone.size() - 1]
            });
        } else {
            this->anchorAtoms.push_back({});
        }
    }
}

std::vector<float> Foldcomp::checkTorsionReconstruction() {
    // Continuize torsion angles
    this->phi = this->phiDisc.continuize(this->phiDiscretized);
    this->psi = this->psiDisc.continuize(this->psiDiscretized);
    this->omega = this->omegaDisc.continuize(this->omegaDiscretized);
    // Append psi, omega, phi to torsion angles
    std::vector<float> output;
    output.reserve(this->phi.size());
    for (size_t i = 0; i < this->phi.size(); i++) {
        output.push_back(this->psi[i]);
        output.push_back(this->omega[i]);
        output.push_back(this->phi[i]);
    }
    return output;
}

int Foldcomp::decompressBackbone(std::vector<AtomCoordinate>& atom) {
    int success;
    atom.clear();
    this->residues.clear();
    this->residueThreeLetter.clear();
    this->backboneTorsionAngles.clear();
    this->backboneBondAngles.clear();
    this->sideChainAnglesPerResidue.clear();
    this->tempFactors.clear();
    this->phi.clear();
    this->psi.clear();
    this->omega.clear();
    this->n_ca_c_angle.clear();
    this->ca_c_n_angle.clear();
    this->c_n_ca_angle.clear();
    this->residues.reserve(this->compressedBackBone.size());
    this->residueThreeLetter.reserve(this->compressedBackBone.size());
    if (this->compressedBackBone.size() > 1) {
        this->backboneTorsionAngles.reserve((this->compressedBackBone.size() - 1) * 3);
    }
    this->backboneBondAngles.reserve(this->compressedBackBone.size() * 3);
    this->sideChainAnglesPerResidue.reserve(this->compressedBackBone.size());
    this->tempFactors.reserve(this->compressedBackBone.size());
    atom.reserve(static_cast<size_t>(this->header.nAtom) + (this->hasOXT ? 1u : 0u));

    // Continuize torsion angles
    this->phi = this->phiDisc.continuize(this->phiDiscretized);
    this->psi = this->psiDisc.continuize(this->psiDiscretized);
    this->omega = this->omegaDisc.continuize(this->omegaDiscretized);

    // Append psi, omega, phi to torsion angles
    for (size_t i = 0; i < this->phi.size(); i++) {
        this->backboneTorsionAngles.push_back(this->psi[i]);
        this->backboneTorsionAngles.push_back(this->omega[i]);
        this->backboneTorsionAngles.push_back(this->phi[i]);
    }

    // Continuize bond angles
    this->n_ca_c_angle = this->n_ca_c_angleDisc.continuize(this->n_ca_c_angleDiscretized);
    this->ca_c_n_angle = this->ca_c_n_angleDisc.continuize(this->ca_c_n_angleDiscretized);
    this->c_n_ca_angle = this->c_n_ca_angleDisc.continuize(this->c_n_ca_angleDiscretized);
    // Append n_ca_c_angle, ca_c_n_angle, c_n_ca_angle to bond angles
    for (size_t i = 0; i < this->n_ca_c_angle.size(); i++) {
        this->backboneBondAngles.push_back(this->ca_c_n_angle[i]);
        this->backboneBondAngles.push_back(this->c_n_ca_angle[i]);
        this->backboneBondAngles.push_back(this->n_ca_c_angle[i]);
    }

    // Get the residue vector
    // TODO: EXTRACT RESIDUE NAMES FROM COMPRESSED BACKBONE
    success = _restoreResidueNames(this->compressedBackBone, this->header, this->residues, this->residueThreeLetter);

    //TODO: FILL IN THIS FUNCTION
    // nerf.reconstruct();
    std::vector<AtomCoordinate> prevForAnchor;
    std::vector<AtomCoordinate> atomByAnchor;
    for (int i = 0; i < this->nAllAnchor - 1; i++) {
        if (i == 0) {
            prevForAnchor = this->prevAtoms;
        }
        // std::vector<int>   sub(&data[100000],&data[101000]);
        // MANUALLY CHECK THAT THE INDICES ARE WITHIN BOUNDS
        size_t compressedSize = this->compressedBackBone.size();
        size_t firstIndex = std::min<size_t>(this->anchorIndices[i], compressedSize);
        size_t lastIndexExclusive = std::min<size_t>(this->anchorIndices[i + 1] + 1, compressedSize);
        std::vector<BackboneChain> subBackbone(
            this->compressedBackBone.begin() + firstIndex,
            this->compressedBackBone.begin() + lastIndexExclusive
        );
        atomByAnchor = reconstructBackboneAtoms(prevForAnchor, subBackbone, this->header);

        // Subset torsion_angles
        size_t torsionSize = this->backboneTorsionAngles.size();
        firstIndex = std::min<size_t>(this->anchorIndices[i] * 3, torsionSize);
        lastIndexExclusive = std::min<size_t>(this->anchorIndices[i + 1] * 3, torsionSize);
        std::vector<float> subTorsionAngles(
            this->backboneTorsionAngles.begin() + firstIndex,
            this->backboneTorsionAngles.begin() + lastIndexExclusive
        );

        success = reconstructBackboneReverse(
            atomByAnchor, this->anchorCoordinates[i], subTorsionAngles, this->nerf
        );
        if (i != this->nAllAnchor - 2) {
            atom.insert(atom.end(), atomByAnchor.begin(), atomByAnchor.end() - 3);
        } else {
            atom.insert(atom.end(), atomByAnchor.begin(), atomByAnchor.end());
        }
        // Update prevForAnchor - last 3 atoms of atomByAnchor
        prevForAnchor = std::vector<AtomCoordinate>(
            atomByAnchor.end() - 3, atomByAnchor.end()
        );
    }

    return success;
}

int Foldcomp::decompress(std::vector<AtomCoordinate>& atom) {
    // 2022-11-15 11:47:49 - Removed defining new vectors for TAs & BAs
    int success = this->decompressBackbone(atom);
    if (success != 0) {
        return success;
    }

    // Reconstruct sidechain
    std::vector<std::pair<size_t, size_t>> residueRanges = splitResidueRanges(atom);
    std::string currResidue = getThreeLetterCode(this->header.firstResidue);

    success = this->_continuizeSideChainTorsionAngles(
        this->sideChainAnglesDiscretized, this->sideChainAnglesPerResidue
    );

    // Prepare tempFactor
    std::vector<float> tempFactors = this->tempFactorsDisc.continuize(this->tempFactorsDiscretized);
    this->tempFactors = tempFactors;

    std::vector<AtomCoordinate> rebuiltAtoms;
    rebuiltAtoms.reserve(static_cast<size_t>(this->header.nAtom) + (this->hasOXT ? 1u : 0u));
    for (size_t i = 0; i < residueRanges.size(); i++) {
        const auto& residueRange = residueRanges[i];
        tcb::span<const AtomCoordinate> backboneResidue(
            atom.data() + residueRange.first, residueRange.second - residueRange.first
        );
        if (i != 0) {
            currResidue = backboneResidue[0].residue;
        }
        const auto aaIt = AAS.find(currResidue);
        std::vector<AtomCoordinate> fullResidue;
        if (aaIt == AAS.end()) {
            fullResidue.assign(backboneResidue.begin(), backboneResidue.end());
        } else {
            fullResidue = nerf.reconstructAminoAcid(
                backboneResidue, this->sideChainAnglesPerResidue[i], aaIt->second
            );
            if (this->useAltAtomOrder) {
                _reorderAtoms(fullResidue, aaIt->second);
            }
        }
        for (AtomCoordinate& residueAtom : fullResidue) {
            residueAtom.tempFactor = tempFactors[i];
            rebuiltAtoms.push_back(std::move(residueAtom));
        }
    }
    atom = std::move(rebuiltAtoms);
    // Reindex atom index of atom
    if (this->hasOXT) {
        // Set OXT tempFactor
        this->OXT.tempFactor = tempFactors.back();
        atom.push_back(this->OXT);
    }
    setAtomIndexSequentially(atom, this->header.idxAtom);

    return success;
}

float Foldcomp::computeBackboneRmsdForCompressedState() const {
    Foldcomp decoded = *this;
    decoded.anchorCoordinates.clear();
    if (decoded.anchorAtoms.size() < 2) {
        return std::numeric_limits<float>::infinity();
    }
    decoded.anchorCoordinates.reserve(decoded.anchorAtoms.size() - 1);
    for (size_t i = 1; i < decoded.anchorAtoms.size(); ++i) {
        const auto& anchorAtoms = decoded.anchorAtoms[i];
        if (anchorAtoms.size() != 3) {
            return std::numeric_limits<float>::infinity();
        }
        std::vector<std::vector<float>> coords;
        coords.reserve(3);
        for (const auto& atom : anchorAtoms) {
            coords.push_back({
                atom.coordinate.x,
                atom.coordinate.y,
                atom.coordinate.z
            });
        }
        decoded.anchorCoordinates.push_back(std::move(coords));
    }

    std::vector<AtomCoordinate> decodedBackbone;
    if (decoded.decompressBackbone(decodedBackbone) != 0) {
        return std::numeric_limits<float>::infinity();
    }
    if (decodedBackbone.size() != this->backbone.size() || decodedBackbone.empty()) {
        return std::numeric_limits<float>::infinity();
    }
    double sumSquaredDistance = 0.0;
    for (size_t i = 0; i < decodedBackbone.size(); ++i) {
        const float dx = this->backbone[i].coordinate.x - decodedBackbone[i].coordinate.x;
        const float dy = this->backbone[i].coordinate.y - decodedBackbone[i].coordinate.y;
        const float dz = this->backbone[i].coordinate.z - decodedBackbone[i].coordinate.z;
        sumSquaredDistance += static_cast<double>(dx) * dx +
                              static_cast<double>(dy) * dy +
                              static_cast<double>(dz) * dz;
    }
    return static_cast<float>(std::sqrt(sumSquaredDistance / decodedBackbone.size()));
}

bool Foldcomp::exceedsBackboneRmsdThreshold(float maxBackboneRmsd) const {
    if (!std::isfinite(maxBackboneRmsd)) {
        return false;
    }
    const double backboneRmsd = computeBackboneRmsdForCompressedState();
    return backboneRmsd > static_cast<double>(maxBackboneRmsd);
}

int Foldcomp::read(const char* data, size_t size) {
    auto readBytes = [&](void* dst, size_t count, size_t& offset) -> bool {
        if (offset + count > size) {
            return false;
        }
        memcpy(dst, data + offset, count);
        offset += count;
        return true;
    };

    this->anchorCoordinates.clear();
    this->lastAtomCoordinates.clear();
    this->compressedBackBone.clear();
    this->sideChainAnglesDiscretized.clear();
    this->tempFactorsDiscretized.clear();

    size_t offset = 0;
    char mNum[MAGICNUMBER_LENGTH];
    if (!readBytes(mNum, MAGICNUMBER_LENGTH, offset)) {
        return -1;
    }
    for (int i = 0; i < MAGICNUMBER_LENGTH; i++) {
        if (mNum[i] != MAGICNUMBER[i]) {
            return -1;
        }
    }

    if (!readBytes(&this->header, sizeof(this->header), offset)) {
        return -1;
    }
    this->read_header(this->header);
    this->anchorIndices.resize(this->nAllAnchor);
    if (!this->anchorIndices.empty() &&
        !readBytes(this->anchorIndices.data(), sizeof(int) * this->nAllAnchor, offset)) {
        return -1;
    }
    this->strTitle.assign(this->header.lenTitle, '\0');
    if (this->header.lenTitle > 0 &&
        !readBytes(this->strTitle.data(), sizeof(char) * this->header.lenTitle, offset)) {
        return -1;
    }

    float prevAtomCoords[9];
    if (!readBytes(prevAtomCoords, sizeof(prevAtomCoords), offset)) {
        return -1;
    }

    if (this->header.nAnchor > 2) {
        float innerAnchorCoords[3];
        for (int i = 0; i < (this->header.nAnchor - 2); i++) {
            std::vector<std::vector<float>> innerAnchorCoord;
            innerAnchorCoord.reserve(3);
            for (int j = 0; j < 3; j++) {
                if (!readBytes(innerAnchorCoords, sizeof(innerAnchorCoords), offset)) {
                    return -1;
                }
                innerAnchorCoord.push_back({
                    innerAnchorCoords[0],
                    innerAnchorCoords[1],
                    innerAnchorCoords[2]
                });
            }
            this->anchorCoordinates.push_back(innerAnchorCoord);
        }
    }

    float lastAtomCoords[9];
    if (!readBytes(lastAtomCoords, sizeof(lastAtomCoords), offset)) {
        return -1;
    }
    for (int i = 0; i < 3; i++) {
        this->lastAtomCoordinates.push_back({ lastAtomCoords[i*3], lastAtomCoords[i*3 + 1], lastAtomCoords[i*3 + 2] });
    }
    this->anchorCoordinates.push_back(this->lastAtomCoordinates);

    if (!readBytes(&this->hasOXT, sizeof(char), offset)) {
        return -1;
    }
    float oxtCoords[3];
    if (!readBytes(oxtCoords, sizeof(oxtCoords), offset)) {
        return -1;
    }
    this->OXT_coords = { oxtCoords[0], oxtCoords[1], oxtCoords[2] };
    int lastResidueIndex = static_cast<int>(this->header.idxResidue + this->header.nResidue - 1);
    int lastAtomIndex = static_cast<int>(this->header.idxAtom + this->header.nAtom - 1);
    this->OXT = AtomCoordinate(
        "OXT", getThreeLetterCode(this->header.lastResidue), getChainName(this->header),
        lastAtomIndex, lastResidueIndex, this->OXT_coords
    );

    this->compressedBackBone.resize(this->header.nResidue);
    char buffer[8];
    for (int i = 0; i < this->header.nResidue; i++) {
        if (!readBytes(buffer, sizeof(buffer), offset)) {
            return -1;
        }
        this->compressedBackBone[i] = convertBytesToBackboneChain(buffer);
    }

    std::vector<unsigned char> encodedSideChain(this->header.nSideChainTorsion);
    if (this->header.nSideChainTorsion > 0 &&
        !readBytes(encodedSideChain.data(), this->header.nSideChainTorsion, offset)) {
        return -1;
    }
    this->sideChainAnglesDiscretized.reserve(this->header.nSideChainTorsion);
    for (size_t i = 0; i < this->header.nSideChainTorsion; i++) {
        this->sideChainAnglesDiscretized.push_back(static_cast<unsigned int>(encodedSideChain[i]));
    }

    if (!readBytes(&this->tempFactorsDisc.min, sizeof(float), offset) ||
        !readBytes(&this->tempFactorsDisc.cont_f, sizeof(float), offset)) {
        return -1;
    }
    std::vector<unsigned char> encodedTempFactors(this->header.nResidue);
    if (this->header.nResidue > 0 &&
        !readBytes(encodedTempFactors.data(), this->header.nResidue, offset)) {
        return -1;
    }
    this->tempFactorsDiscretized.reserve(this->header.nResidue);
    for (int i = 0; i < this->header.nResidue; i++) {
        this->tempFactorsDiscretized.push_back(static_cast<unsigned int>(encodedTempFactors[i]));
    }

    int success = _restoreAtomCoordinate(prevAtomCoords);
    if (success != 0) {
        return -2;
    }
    for (int i = 0; i < 6; i++) {
        success = _restoreDiscretizer(i);
        if (success != 0) {
            return success;
        }
    }
    return success;
}

int Foldcomp::writeString(std::string& output) {
    int flag = 0;
    output.clear();
    size_t totalSize = MAGICNUMBER_LENGTH +
                       sizeof(CompressedFileHeader) +
                       this->anchorIndices.size() * sizeof(int) +
                       this->strTitle.size() +
                       this->anchorAtoms.size() * 9 * sizeof(float) +
                       sizeof(char) +
                       3 * sizeof(float) +
                       this->compressedBackBone.size() * 8 +
                       this->sideChainAnglesDiscretized.size() +
                       2 * sizeof(float) +
                       this->tempFactorsDiscretized.size();
    output.reserve(totalSize);

    output.append(MAGICNUMBER, MAGICNUMBER_LENGTH);
    appendBytes(output, this->header);
    appendBytes(output, this->anchorIndices);
    output.append(this->strTitle);

    for (const auto& anchors : this->anchorAtoms) {
        for (int i = 0; i < 3; i++) {
            float3d coord;
            if ((size_t)i < anchors.size()) {
                coord = anchors[i].coordinate;
            } else if (!anchors.empty()) {
                coord = anchors.back().coordinate;
            } else {
                coord = {0.0f, 0.0f, 0.0f};
            }
            appendBytes(output, coord.x);
            appendBytes(output, coord.y);
            appendBytes(output, coord.z);
        }
    }

    output.push_back(this->hasOXT);
    appendBytes(output, this->OXT_coords.x);
    appendBytes(output, this->OXT_coords.y);
    appendBytes(output, this->OXT_coords.z);

    char buffer[8];
    for (const auto& residue : this->compressedBackBone) {
        BackboneChain tmp = residue;
        flag = convertBackboneChainToBytes(tmp, buffer);
        output.append(buffer, sizeof(buffer));
    }

    appendByteVector(output, this->sideChainAnglesDiscretized);
    appendBytes(output, this->tempFactorsDisc.min);
    appendBytes(output, this->tempFactorsDisc.cont_f);
    appendByteVector(output, this->tempFactorsDiscretized);
    return flag;
}

int Foldcomp::write(std::string filename) {
    std::string encoded;
    int flag = writeString(encoded);
    if (flag != 0) {
        return flag;
    }
    FILE* outfile = fopen(filename.c_str(), "wb");
    if (outfile == nullptr) {
        return -1;
    }
    size_t written = fwrite(encoded.data(), 1, encoded.size(), outfile);
    fclose(outfile);
    return written == encoded.size() ? 0 : -1;
}

#ifdef FOLDCOMP_EXECUTABLE
// 2022-08-29 15:42:29 TAR format support
int Foldcomp::writeTar(mtar_t& tar, std::string filename, size_t size) {
    int flag = 0;
    mtar_write_file_header(&tar, filename.c_str(), size);
    // Magic number
    mtar_write_data(&tar, MAGICNUMBER, MAGICNUMBER_LENGTH);
    // Write header
    mtar_write_data(&tar, &this->header, sizeof(CompressedFileHeader));
    // Write anchorIndices
    for (size_t i = 0; i < this->anchorIndices.size(); i++) {
        mtar_write_data(&tar, &this->anchorIndices[i], sizeof(int));
    }
    // Write title
    mtar_write_data(&tar, this->strTitle.c_str(), this->strTitle.length());
    // Write anchor atoms
    for (const auto& anchors : this->anchorAtoms) {
        for (int i = 0; i < 3; i++) {
            float3d coord;
            if ((size_t)i < anchors.size()) {
                coord = anchors[i].coordinate;
            } else if (!anchors.empty()) {
                coord = anchors.back().coordinate;
            } else {
                coord = {0.0f, 0.0f, 0.0f};
            }
            mtar_write_data(&tar, &coord.x, sizeof(float));
            mtar_write_data(&tar, &coord.y, sizeof(float));
            mtar_write_data(&tar, &coord.z, sizeof(float));
        }
    }
    // Write hasOXT
    mtar_write_data(&tar, &this->hasOXT, sizeof(char));
    // Write OXT_coords
    mtar_write_data(&tar, &this->OXT_coords.x, sizeof(float));
    mtar_write_data(&tar, &this->OXT_coords.y, sizeof(float));
    mtar_write_data(&tar, &this->OXT_coords.z, sizeof(float));
    // Write sideChainDisc
    // mtar_write_data(&tar, &this->sideChainDisc, sizeof(SideChainDiscretizers));
    // Write the compressed backbone
    char* buffer = new char[8];
    for (size_t i = 0; i < this->compressedBackBone.size(); i++) {
        flag = convertBackboneChainToBytes(this->compressedBackBone[i], buffer);
        mtar_write_data(&tar, buffer, 8);
    }
    delete[] buffer;
    // Write side chain torsion angles
    int encodedSideChainSize = this->nSideChainTorsion;
    if (encodedSideChainSize % 2 == 1) {
        encodedSideChainSize++;
    }
    encodedSideChainSize /= 2;

    unsigned char* charSideChainTorsion = new unsigned char[this->nSideChainTorsion];
    // Get array of unsigned int from sideChainAnglesDiscretized and convert to char array
    for (int i = 0; i < this->nSideChainTorsion; i++) {
        // convert unsigned int to char
        charSideChainTorsion[i] = (unsigned char)this->sideChainAnglesDiscretized[i];
    }
    mtar_write_data(&tar, charSideChainTorsion, this->sideChainAnglesDiscretized.size());
    delete[] charSideChainTorsion;
    // Write temperature factors
    // Disc
    mtar_write_data(&tar, &this->tempFactorsDisc.min, sizeof(float));
    mtar_write_data(&tar, &this->tempFactorsDisc.cont_f, sizeof(float));
    // Convert unsigned int to char array
    unsigned char* charTempFactors = new unsigned char[this->header.nResidue];
    // Get array of unsigned int from tempFactorsDiscretized and convert to char array
    for (int i = 0; i < this->header.nResidue; i++) {
        // convert unsigned int to char
        charTempFactors[i] = (unsigned char)this->tempFactorsDiscretized[i];
    }
    mtar_write_data(&tar, charTempFactors, this->tempFactorsDiscretized.size());
    delete[] charTempFactors;
    return flag;
}
#endif

size_t Foldcomp::getSize() {
    // Calculate the size of the compressed format
    size_t size = 0;
    // Magic number
    size += MAGICNUMBER_LENGTH;
    // Header
    size += sizeof(CompressedFileHeader);
    // Anchor indices
    size += sizeof(int) * this->anchorIndices.size();
    // Title
    size += this->strTitle.length();
    // Anchor atoms
    size += sizeof(float) * 3 * 3 * this->anchorAtoms.size();
    // OXT
    size += sizeof(char);
    size += sizeof(float) * 3;
    // Backbone
    size += sizeof(char) * 8 * this->compressedBackBone.size();
    // Side chain torsion angles
    size += sizeof(unsigned char) * this->sideChainAnglesDiscretized.size();
    // Temperature factors
    size += sizeof(float) * 2;
    size += sizeof(unsigned char) * this->tempFactorsDiscretized.size();
    return size;
}


// Functions to extract temperature factors only
int Foldcomp::continuizeTempFactors() {
    this->tempFactors = this->tempFactorsDisc.continuize(this->tempFactorsDiscretized);
    return 0;
}

int Foldcomp::writeTorsionAngles(std::string filename) {
    FILE* outfile = fopen(filename.c_str(), "wb");
    if (outfile == nullptr) {
        return -1;
    }
    fputs("index,phi,psi,omega\n", outfile);
    for (size_t i = 0; i < this->phi.size(); i++) {
        fprintf(outfile, "%zu,%f,%f,%f\n", i, this->phi[i], this->psi[i], this->omega[i]);
    }
    fclose(outfile);
    return 0;
}

/**
 * @brief Extract information from the compressed file and write to a FASTA-like file
 *
 * @param filename
 * @param type 0: plddt, 1: sequence
 * @return int
 */
int Foldcomp::extract(std::string& data, int type, int digits) {
    int flag = 0;
    if (type == 0) {
        // Extract temperature factors
        this->continuizeTempFactors();
        if (digits < 1) {
            digits = 1;
        } else if (digits > 4) {
            digits = 4;
        }
        // Reserve string size
        int reserving;
        if (digits == 1) {
            reserving = this->tempFactors.size();
        } else if (digits == 2) {
            // 2 digits, comma separated
            reserving = this->tempFactors.size() * 3;
        } else if (digits == 3) {
            // 3 digits, one decimal place, comma separated
            reserving = this->tempFactors.size() * 5;
        } else {
            // 4 digits, two decimal places, comma separated
            reserving = this->tempFactors.size() * 6;
        }
        data.reserve(reserving);
        // Check if the max value is 0.0-1.0 or 1-100
        bool isZeroToOne = false;
        float maxval = this->tempFactorsDisc.cont_f * (pow(2, NUM_BITS_TEMP) - 1) + this->tempFactorsDisc.min;
        if (maxval <= 1.0f && digits <= 2) {
            isZeroToOne = true;
        }

        for (size_t i = 0; i < this->tempFactors.size(); i++) {
            float clamped;
            char digit1, digit2;

            if (isZeroToOne) {
                clamped = (std::clamp(this->tempFactors[i], 0.0f, 1.0f));
                digit1 = (char)((int)(clamped * 10.0f) % 10) + '0';
                digit2 = (char)((int)(clamped * 100.0f) % 10) + '0';
            } else {
                clamped = (std::clamp(this->tempFactors[i], 0.0f, 100.0f));
                digit1 = (char)(clamped / 10.0f) + '0';
                digit2 = (char)((int)clamped % 10) + '0';
            }
            // Append digits to data
            data.append(1, digit1);
            if (digits > 1) {
                data.append(1, digit2);
            }

            if (digits >= 3) {
                char digit3 = (char)((int)(clamped * 10.0f) % 10) + '0';
                data.append(1, '.');
                data.append(1, digit3);
            }

            if (digits == 4) {
                char digit4 = (char)((int)(clamped * 100.0f) % 10) + '0';
                data.append(1, digit4);
            }

            if (digits > 1 && i != this->tempFactors.size() - 1) {
                // Add comma separator
                data.append(1, ',');
            }
        }
    } else if (type == 1) {
        // Extract sequence
        data.reserve(this->header.nResidue);
        for (int i = 0; i < this->header.nResidue; i++) {
            char res = convertIntToOneLetterCode(this->compressedBackBone[i].residue);
            data.append(1, res);
        }
    }
    return flag;
}


//
CompressedFileHeader Foldcomp::get_header() {
    CompressedFileHeader header = {};
    // counts
    header.nResidue = this->nResidue;
    header.nAtom = this->nAtom;
    header.idxResidue = this->idxResidue;
    header.idxAtom = this->idxAtom;
    header.nAnchor = this->nAllAnchor;
    header.version = 0;
    header.flags = 0;
    header.nSideChainTorsion = this->nSideChainTorsion;
    header.firstResidue = this->firstResidue;
    header.lastResidue = this->lastResidue;
    header.lenTitle = this->lenTitle;
    std::string chainName;
    if (!this->backbone.empty()) {
        chainName = this->backbone[0].chain;
    } else if (this->chain != '\0') {
        chainName.assign(1, this->chain);
    }
    setChainName(header, chainName);
    // discretizer parameters
    header.mins[0] = this->phiDisc.min;
    header.mins[1] = this->psiDisc.min;
    header.mins[2] = this->omegaDisc.min;
    header.mins[3] = this->n_ca_c_angleDisc.min;
    header.mins[4] = this->ca_c_n_angleDisc.min;
    header.mins[5] = this->c_n_ca_angleDisc.min;
    header.cont_fs[0] = this->phiDisc.cont_f;
    header.cont_fs[1] = this->psiDisc.cont_f;
    header.cont_fs[2] = this->omegaDisc.cont_f;
    header.cont_fs[3] = this->n_ca_c_angleDisc.cont_f;
    header.cont_fs[4] = this->ca_c_n_angleDisc.cont_f;
    header.cont_fs[5] = this->c_n_ca_angleDisc.cont_f;
    return header;
}

int Foldcomp::read_header(CompressedFileHeader& header) {
    this->nResidue = header.nResidue;
    this->nAtom = header.nAtom;
    this->idxResidue = header.idxResidue;
    this->idxAtom = header.idxAtom;
    this->nAllAnchor = header.nAnchor;
    this->nSideChainTorsion = header.nSideChainTorsion;
    this->firstResidue = header.firstResidue;
    this->lastResidue = header.lastResidue;
    this->lenTitle = header.lenTitle;
    //
    this->chain = header.chain;
    // discretizer parameters
    this->phiDisc.min = header.mins[0];
    this->psiDisc.min = header.mins[1];
    this->omegaDisc.min = header.mins[2];
    this->n_ca_c_angleDisc.min = header.mins[3];
    this->ca_c_n_angleDisc.min = header.mins[4];
    this->c_n_ca_angleDisc.min = header.mins[5];
    this->phiDisc.cont_f = header.cont_fs[0];
    this->psiDisc.cont_f = header.cont_fs[1];
    this->omegaDisc.cont_f = header.cont_fs[2];
    this->n_ca_c_angleDisc.cont_f = header.cont_fs[3];
    this->ca_c_n_angleDisc.cont_f = header.cont_fs[4];
    this->c_n_ca_angleDisc.cont_f = header.cont_fs[5];
    return 0;
}

void Foldcomp::print(int length) {
    // Print the header
    std::cout << "[Header]" << std::endl;
    std::cout << "nResidue: " << this->header.nResidue << std::endl;
    std::cout << "nAtom: " << this->header.nAtom << std::endl;
    std::cout << "idxResidue: " << this->header.idxResidue << std::endl;
    std::cout << "idxAtom: " << this->header.idxAtom << std::endl;
    std::cout << "nSideChainTorsion: " << this->header.nSideChainTorsion << std::endl;
    std::cout << "mins: " << std::endl;
    for (int i = 0; i < 6; i++) {
        std::cout << this->header.mins[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "cont_fs: " << std::endl;
    for (int i = 0; i < 6; i++) {
        std::cout << this->header.cont_fs[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "--------------------" << std::endl;

    // Print the prevAtoms
    std::cout << "[PrevAtoms]" << std::endl;
    for (int i = 0; i < 3; i++) {
        std::cout << "Atom " << this->prevAtoms[i].atom << ": " << std::endl;
        std::cout << "x: " << this->prevAtoms[i].coordinate.x << std::endl;
        std::cout << "y: " << this->prevAtoms[i].coordinate.y << std::endl;
        std::cout << "z: " << this->prevAtoms[i].coordinate.z << std::endl;
    }
    std::cout << "--------------------" << std::endl;

    // Print the first element of compressedBackBone
    std::cout << "[CompressedBackbone]" << std::endl;
    for (int i = 0; i < length; i++) {
        std::cout << "Residue: " << this->compressedBackBone[i].residue << std::endl;
        std::cout << "phi-disc: " << this->compressedBackBone[i].phi << " / ";
        std::cout << this->phiDisc.continuize(this->compressedBackBone[i].phi) << std::endl;
        std::cout << "psi-disc: " << this->compressedBackBone[i].psi << " / ";
        std::cout << this->psiDisc.continuize(this->compressedBackBone[i].psi) << std::endl;
        std::cout << "omega-disc: " << this->compressedBackBone[i].omega << " / ";
        std::cout << this->omegaDisc.continuize(this->compressedBackBone[i].omega) << std::endl;
        std::cout << "n_ca_c_angle-disc: " << this->compressedBackBone[i].n_ca_c_angle << " / ";
        std::cout << this->n_ca_c_angleDisc.continuize(this->compressedBackBone[i].n_ca_c_angle) << std::endl;
        std::cout << "ca_c_n_angle-disc: " << this->compressedBackBone[i].ca_c_n_angle << " / ";
        std::cout << this->ca_c_n_angleDisc.continuize(this->compressedBackBone[i].ca_c_n_angle) << std::endl;
        std::cout << "c_n_ca_angle-disc: " << this->compressedBackBone[i].c_n_ca_angle << " / ";
        std::cout << this->c_n_ca_angleDisc.continuize(this->compressedBackBone[i].c_n_ca_angle) << std::endl;
    }
}

/*
void Foldcomp::printSideChainTorsion(std::string filename) {
    FILE* outfile = fopen(filename.c_str(), "wb");
    //
    outfile << "ResidueInd,Residue,Type,Key,RawVal,DiscVal,DiscMin,DiscContF,ReconVal,Diff\n";
    int movingIndex = 0;
    std::vector< std::vector<AtomCoordinate> > atomByResidue = splitAtomByResidue(this->rawAtoms);
    std::map<std::string, float> currBondAngle;
    std::map<std::string, float> currBondLength;
    std::map<std::string, float> currTorsionAngle;

    for (int i = 0; i < this->nResidue; i++) {
        currBondAngle = calculateBondAngles(atomByResidue[i], this->AAS.at(this->residueThreeLetter[i]));
        currBondLength = calculateBondLengths(atomByResidue[i], this->AAS.at(this->residueThreeLetter[i]));
        for (const auto& bl: currBondLength) {
            outfile << i << "," << this->residueThreeLetter[i] << ",BondLength,";
            outfile << bl.first << "," << bl.second << ",NA,NA,NA,";
            outfile << this->AAS.at(this->residueThreeLetter[i]).bondLengths.at(bl.first) << ",";
            outfile << bl.second - this->AAS.at(this->residueThreeLetter[i]).bondLengths.at(bl.first) << "\n";
        }
        for (const auto& ba : currBondAngle) {
            outfile << i << "," << this->residueThreeLetter[i] << ",BondAngle,";
            outfile << ba.first << "," << ba.second << ",NA,NA,NA," << this->AAS.at(this->residueThreeLetter[i]).bondAngles.at(ba.first) << ",";
            outfile << ba.second - this->AAS.at(this->residueThreeLetter[i]).bondAngles.at(ba.first) << "\n";
        }
        for (size_t j = 0; j < this->sideChainAnglesPerResidue[i].size(); j++) {
            outfile << i << "," << this->residueThreeLetter[i] << ",TorsionAngle,";
            outfile << j << "," << this->sideChainAnglesPerResidue[i][j] << ",";
            outfile << this->sideChainAnglesDiscretized[movingIndex] << ",";
            outfile << this->sideChainDiscMap[this->residueThreeLetter[i]][j].min << ",";
            outfile << this->sideChainDiscMap[this->residueThreeLetter[i]][j].cont_f << ",";
            outfile << this->sideChainDiscMap[this->residueThreeLetter[i]][j].continuize(this->sideChainAnglesDiscretized[movingIndex]) << ",";
            outfile << this->sideChainAnglesPerResidue[i][j] - this->sideChainDiscMap[this->residueThreeLetter[i]][j].continuize(this->sideChainAnglesDiscretized[movingIndex]) << "\n";
            movingIndex++;
        }
    }
    outfile.close();
}
*/

/**
 * @brief Checks the input file read is valid or not.
 *        This method is expected to be called after the input file is read.
 * @return int Error code. 0 if no error.
 */
ValidityError Foldcomp::checkValidity() {
    // Check size
    bool hasCorrectNumResidue = (this->header.nResidue == this->compressedBackBone.size());
    bool hasCorrectNumSideChain = (this->header.nSideChainTorsion == this->sideChainAnglesDiscretized.size());
    bool hasCorrectNumTempfactor = (this->header.nResidue == this->tempFactorsDiscretized.size());
    // Check vectors are not empty
    // For backbone, just check the torsion angles
    bool emptyBackbone = std::all_of(
        this->compressedBackBone.begin(),
        this->compressedBackBone.end(),
        [](BackboneChain bb) {
            return bb.phi == 0 && bb.psi == 0 && bb.omega == 0;
        }
    );
    bool emptySideChain = std::all_of(
        this->sideChainAnglesDiscretized.begin(),
        this->sideChainAnglesDiscretized.end(),
        [](unsigned int i){return i == 0;}
    );
    bool emptyTempFactor = std::all_of(
        this->tempFactorsDiscretized.begin(),
        this->tempFactorsDiscretized.end(),
        [](unsigned int i){return i == 0;}
    );
    // Return
    if (!hasCorrectNumResidue) {
        return E_BACKBONE_COUNT_MISMATCH;
    } else if (!hasCorrectNumSideChain) {
        return E_SIDECHAIN_COUNT_MISMATCH;
    } else if (!hasCorrectNumTempfactor) {
        return E_TEMP_FACTOR_COUNT_MISMATCH;
    } else if (emptyBackbone) {
        return E_EMPTY_BACKBONE_ANGLE;
    } else if (emptySideChain) {
        return E_EMPTY_SIDECHAIN_ANGLE;
    } else if (emptyTempFactor) {
        return E_EMPTY_TEMP_FACTOR;
    } else {
        return SUCCESS;
    }
}

void printValidityError(ValidityError err, std::string& filename) {
    // Print error message to stderr with filename
    switch (err) {
        case E_BACKBONE_COUNT_MISMATCH:
            std::clog << "[Error] Number of backbone angles does not match header: " << filename << std::endl;
            break;
        case E_SIDECHAIN_COUNT_MISMATCH:
            std::clog << "[Error] Number of sidechain angles does not match header: " << filename << std::endl;
            break;
        case E_TEMP_FACTOR_COUNT_MISMATCH:
            std::clog << "[Error] Number of temperature factors does not match header: " << filename << std::endl;
            break;
        case E_EMPTY_BACKBONE_ANGLE:
            std::clog << "[Error] All backbone angles are empty: " << filename << std::endl;
            break;
        case E_EMPTY_SIDECHAIN_ANGLE:
            std::clog << "[Error] All sidechain angles are empty: " << filename << std::endl;
            break;
        case E_EMPTY_TEMP_FACTOR:
            std::clog << "[Error] All temperature factors are empty: " << filename << std::endl;
            break;
        case SUCCESS:
            break;
        default:
            std::clog << "[Error] Unknown error: " << filename << std::endl;
            break;
    }
}

void _reorderAtoms(std::vector<AtomCoordinate>& atoms, const AminoAcid& aa) {
    std::vector<AtomCoordinate> newAtoms = atoms;
    for (size_t i = 0; i < atoms.size(); i++) {
        if (atoms[i].atom == aa.altAtoms[i]) {
            continue;
        } else {
            for (size_t j = 0; j < aa.altAtoms.size(); j++) {
                if (atoms[i].atom == aa.altAtoms[j]) {
                    newAtoms[j] = atoms[i];
                }
            }
        }
    }
    atoms = newAtoms;
}


// Sidechain

char* encodeSideChainTorsionVector(std::vector<unsigned int> vector) {
    //
    size_t size = vector.size();
    size_t newSize = size;
    if (size % 2 == 1) {
        newSize++;
    }
    newSize /= 2;
    char* output = new char[size];
    char temp;
    for (size_t i = 0; i < newSize; i++) {
        // First 4 bits
        temp = vector[i * 2] & 0x0F;
        temp = temp << 4;
        if (i * 2 + 1 < vector.size()) {
            temp |= (vector[i * 2 + 1] & 0x0F);
        } else {
            temp |= 0x0F;
        }
        output[i] = temp;
    }
    return output;
}

int decodeSideChainTorsionVector(char* input, int nTorsion, std::vector<unsigned int>& vector) {
    int size = nTorsion;
    if (size % 2 == 1) {
        size++;
    }
    size /= 2;
    unsigned int temp, first, second;
    if (!vector.empty()) {
        vector.clear();
    }
    for (int i = 0; i < size; i++) {
        temp = input[i];
        first = temp >> 4;
        second = temp & 0x0F;
        vector.push_back(first);
        if (i * 2 + 1 < nTorsion) {
            vector.push_back(second);
        }
    }
    return size;
}

unsigned char* encodeDiscretizedTempFactors(std::vector<unsigned int> vector) {
    unsigned char* output = new unsigned char[vector.size()];
    for (size_t i = 0; i < vector.size(); i++) {
        output[i] = (unsigned char)vector[i];
    }
    return output;
}

int decodeDiscretizedTempFactors(unsigned char* input, int size, std::vector<unsigned int>& vector) {
    int out = 0;
    unsigned int temp;
    // If vector is not empty, clear it
    if (!vector.empty()) {
        vector.clear();
    }
    for (int i = 0; i < size; i++) {
        temp = (unsigned int)input[i];
        vector.push_back(temp);
    }
    return out;
}

std::map<std::string, std::vector<Discretizer> > initializeSideChainDiscMap() {
    // Initialize the map
    std::map<std::string, std::vector<Discretizer> > discMap;
    // Get Amino acids list
    std::vector<std::string> aaNames = getAminoAcidList();
    int numTorsion = 0;
    // We can access to the specific angle discretizer with AA name and torsion index
    // ex) discmap["ALA"][0]
    for (size_t i = 0; i < aaNames.size(); i++) {
        numTorsion = getSideChainTorsionNum(aaNames[i]);
        discMap[aaNames[i]] = std::vector<Discretizer>(numTorsion);
    }
    return discMap;
}

float* getMinPointerFromSideChainDiscretizers(
    std::string& residue, SideChainDiscretizers& scDiscretizers
) {
    if (residue == "ALA") {
        return scDiscretizers.ala_min;
    } else if (residue == "ARG") {
        return scDiscretizers.arg_min;
    } else if (residue == "ASN") {
        return scDiscretizers.asn_min;
    } else if (residue == "ASP") {
        return scDiscretizers.asp_min;
    } else if (residue == "CYS") {
        return scDiscretizers.cys_min;
    } else if (residue == "GLN") {
        return scDiscretizers.gln_min;
    } else if (residue == "GLU") {
        return scDiscretizers.glu_min;
    } else if (residue == "GLY") {
        return scDiscretizers.gly_min;
    } else if (residue == "HIS") {
        return scDiscretizers.his_min;
    } else if (residue == "ILE") {
        return scDiscretizers.ile_min;
    } else if (residue == "LEU") {
        return scDiscretizers.leu_min;
    } else if (residue == "LYS") {
        return scDiscretizers.lys_min;
    } else if (residue == "MET") {
        return scDiscretizers.met_min;
    } else if (residue == "PHE") {
        return scDiscretizers.phe_min;
    } else if (residue == "PRO") {
        return scDiscretizers.pro_min;
    } else if (residue == "SER") {
        return scDiscretizers.ser_min;
    } else if (residue == "THR") {
        return scDiscretizers.thr_min;
    } else if (residue == "TRP") {
        return scDiscretizers.trp_min;
    } else if (residue == "TYR") {
        return scDiscretizers.tyr_min;
    } else if (residue == "VAL") {
        return scDiscretizers.val_min;
    } else {
        return NULL;
    }
}

float* getContFFromSideChainDiscretizers(
    std::string& residue, SideChainDiscretizers& scDiscretizers
) {
    if (residue == "ALA") {
        return scDiscretizers.ala_cont_fs;
    } else if (residue == "ARG") {
        return scDiscretizers.arg_cont_fs;
    } else if (residue == "ASN") {
        return scDiscretizers.asn_cont_fs;
    } else if (residue == "ASP") {
        return scDiscretizers.asp_cont_fs;
    } else if (residue == "CYS") {
        return scDiscretizers.cys_cont_fs;
    } else if (residue == "GLN") {
        return scDiscretizers.gln_cont_fs;
    } else if (residue == "GLU") {
        return scDiscretizers.glu_cont_fs;
    } else if (residue == "GLY") {
        return scDiscretizers.gly_cont_fs;
    } else if (residue == "HIS") {
        return scDiscretizers.his_cont_fs;
    } else if (residue == "ILE") {
        return scDiscretizers.ile_cont_fs;
    } else if (residue == "LEU") {
        return scDiscretizers.leu_cont_fs;
    } else if (residue == "LYS") {
        return scDiscretizers.lys_cont_fs;
    } else if (residue == "MET") {
        return scDiscretizers.met_cont_fs;
    } else if (residue == "PHE") {
        return scDiscretizers.phe_cont_fs;
    } else if (residue == "PRO") {
        return scDiscretizers.pro_cont_fs;
    } else if (residue == "SER") {
        return scDiscretizers.ser_cont_fs;
    } else if (residue == "THR") {
        return scDiscretizers.thr_cont_fs;
    } else if (residue == "TRP") {
        return scDiscretizers.trp_cont_fs;
    } else if (residue == "TYR") {
        return scDiscretizers.tyr_cont_fs;
    } else if (residue == "VAL") {
        return scDiscretizers.val_cont_fs;
    } else {
        return NULL;
    }
}

int getSideChainTorsionNum(std::string residue) {
    int out = 0;
    if (residue == "ALA") {
        out = 2;
    } else if (residue == "ARG") {
        out = 8;
    } else if (residue == "ASN") {
        out = 5;
    } else if (residue == "ASP") {
        out = 5;
    } else if (residue == "CYS") {
        out = 3;
    } else if (residue == "GLN") {
        out = 6;
    } else if (residue == "GLU") {
        out = 6;
    } else if (residue == "GLY") {
        out = 1;
    } else if (residue == "HIS") {
        out = 7;
    } else if (residue == "ILE") {
        out = 5;
    } else if (residue == "LEU") {
        out = 5;
    } else if (residue == "LYS") {
        out = 6;
    } else if (residue == "MET") {
        out = 5;
    } else if (residue == "PHE") {
        out = 8;
    } else if (residue == "PRO") {
        out = 4;
    } else if (residue == "SER") {
        out = 3;
    } else if (residue == "THR") {
        out = 4;
    } else if (residue == "TRP") {
        out = 11;
    } else if (residue == "TYR") {
        out = 9;
    } else if (residue == "VAL") {
        out = 4;
    } else {
        out = 0;
    }
    return out;
}

const std::map<std::string, AminoAcid> Foldcomp::AAS = AminoAcid::AminoAcids();
