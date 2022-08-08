/**
 * File: compressed_torsion.cpp
 * Project: src
 * Created: 2021-02-04 13:31:52
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     This file contains main data structures for torsion angle compression and
 *     functions for handling them.
 * ---
 * Last Modified: 2022-08-08 22:31:50
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */

#include "compressed_torsion.h"

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
    unsigned int r = convertOneLetterCodeToInt(residue);
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
    BackboneChain& bb, CompressedFileHeader& header
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
    std::vector<BackboneChain>& bbv, CompressedFileHeader& header
) {
    std::vector<DecompressedBackboneChain> output;
    for (auto& bb : bbv) {
        output.push_back(decompressBackboneChain(bb, header));
    }
    return output;
}

float _continuize(unsigned int input, float min, float cont_f) {
    float output = min + ((float)input * cont_f);
    return output;
}

/**
 * @brief Reconstruct backbone atoms from compressed backbone info
 *
 * @param prev_atoms std::vector<AtomCoordinates> of 3 previous atoms
 * @param backbone std::vector<BackboneChain>
 * @return std::vector<AtomCoordinate>
 */
std::vector<AtomCoordinate> reconstructBackboneAtoms(
    std::vector<AtomCoordinate>& prevAtoms,
    std::vector<BackboneChain>& backbone,
    CompressedFileHeader& header
) {
    Nerf nerf;
    std::vector<AtomCoordinate> reconstructedAtoms;

    // Save first three atoms
    for (int i = 0; i < 3; i++) {
        reconstructedAtoms.push_back(prevAtoms[i]);
    }
    int total = backbone.size();

    std::vector< std::vector<float> > prevCoords;
    std::vector<DecompressedBackboneChain> deBackbone;
    deBackbone = decompressBackboneChain(backbone, header);
    AtomCoordinate currN, currCA, currC;
    std::vector<float> currNCoord, currCACoord, currCCoord;
    float currBondLength, currBondAngle, currTorsionAngle;
    int currAtomIndex = prevAtoms[2].atom_index + 1;
    int currResidueIndex = prevAtoms[2].residue_index + 1;
    std::string currResidue;
    char currResidueChar;

    // Iterate through backbone
    // Should put N, CA, C in this loop
    for (int i = 0; i < (total - 1); i++) {
        prevAtoms = {reconstructedAtoms[i*3], reconstructedAtoms[i*3+1], reconstructedAtoms[i*3+2]};
        prevCoords = extractCoordinates(prevAtoms);
        // Convert char (deBackbone[i].residue) to string (currResidue)
        currResidue = getThreeLetterCode(deBackbone[i + 1].residue);

        // Place N
        currNCoord = nerf.place_atom(
            prevCoords, C_TO_N_DIST, deBackbone[i].ca_c_n_angle, deBackbone[i].psi
        );
        currN = AtomCoordinate(
            "N", currResidue, prevAtoms[0].chain,
            currAtomIndex, currResidueIndex, currNCoord
        );
        // Place CA
        currAtomIndex++;
        prevCoords = {prevCoords[1], prevCoords[2], currNCoord};
        if (deBackbone[i].residue != 'P') {
            currCACoord = nerf.place_atom(
                prevCoords, N_TO_CA_DIST, deBackbone[i].c_n_ca_angle, deBackbone[i].omega
            );
        } else {
            currCACoord = nerf.place_atom(
                prevCoords, PRO_N_TO_CA_DIST, deBackbone[i].c_n_ca_angle, deBackbone[i].omega
            );
        }
        currCA = AtomCoordinate(
            "CA", currResidue, prevAtoms[0].chain,
            currAtomIndex, currResidueIndex, currCACoord
        );
        // Place C
        currAtomIndex++;
        prevCoords = {prevCoords[1], prevCoords[2], currCACoord};
        //
        currCCoord = nerf.place_atom(
            prevCoords, CA_TO_C_DIST, deBackbone[i].n_ca_c_angle, deBackbone[i].phi
        );
        currC = AtomCoordinate(
            "C", currResidue, prevAtoms[0].chain,
            currAtomIndex, currResidueIndex, currCCoord
        );

        // Append to reconstructedAtoms
        reconstructedAtoms.push_back(currN);
        reconstructedAtoms.push_back(currCA);
        reconstructedAtoms.push_back(currC);
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
    std::vector<AtomCoordinate> atomBack = atom;
    // Last atoms
    atomBack[atomBack.size() - 3].coordinate[0] = lastCoords[0][0];
    atomBack[atomBack.size() - 3].coordinate[1] = lastCoords[0][1];
    atomBack[atomBack.size() - 3].coordinate[2] = lastCoords[0][2];
    atomBack[atomBack.size() - 2].coordinate[0] = lastCoords[1][0];
    atomBack[atomBack.size() - 2].coordinate[1] = lastCoords[1][1];
    atomBack[atomBack.size() - 2].coordinate[2] = lastCoords[1][2];
    atomBack[atomBack.size() - 1].coordinate[0] = lastCoords[2][0];
    atomBack[atomBack.size() - 1].coordinate[1] = lastCoords[2][1];
    atomBack[atomBack.size() - 1].coordinate[2] = lastCoords[2][2];

   std::vector<float> bond_angles = nerf.getBondAngles(atom);

    std::vector<AtomCoordinate> atomBackward = nerf.reconstructWithReversed(
        atomBack, torsion_angles, bond_angles
    );

    atomBack = weightedAverage(atom, atomBackward);
    atom = atomBack;
    return 0;
}


int discretizeSideChainTorsionAngles(
    std::vector< std::vector<float> >& torsionPerResidue,
    std::vector<std::string>& residueNames,
    std::map<std::string, AminoAcid>& AAS,
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
    for (auto& sc : sideChainTorsionMap) {
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
    for (int i = 0; i < torsionPerResidue.size(); i++) {
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
    // Declare
    int success = 0;
    std::string currResidue;
    int currResidueTorsionNum = 0;
    int currIndex = 0;
    float currTorsion = 0;

    scDiscretizersMap = initializeSideChainDiscMap();
    success = fillSideChainDiscretizerMap(scDiscretizers, scDiscretizersMap);
    std::vector< std::vector<float> > torsionPerResidue;
    std::vector<float> currTorsionVector;
    FixedAngleDiscretizer currDiscretizer = FixedAngleDiscretizer(pow(2, NUM_BITS_TEMP) - 1);
    // Iterate
    for (int i = 0; i < residueNames.size(); i++) {
        currResidue = residueNames[i];
        currResidueTorsionNum = getSideChainTorsionNum(currResidue);
        currTorsionVector.clear();
        currTorsionVector.resize(currResidueTorsionNum);
        for (int j = 0; j < currResidueTorsionNum; j++) {
            // Get current torsion angle vector
            //currTorsion = scDiscretizersMap[currResidue][j].continuize(torsionDiscretized[currIndex]);
            currTorsion = currDiscretizer.continuize(torsionDiscretized[currIndex]);
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
    for (auto& sc : scDiscretizersMap) {
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

int CompressedResidue::_discretizeSideChainTorsionAngles(
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

int CompressedResidue::_continuizeSideChainTorsionAngles(
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
    int flag = 0;
    flag = convertBackboneChainToBytes(res, byteArray);
    for (int i = 0; i < 8; i++) {
        bits = byteArray[i];
        std::cout << bits << " ";
    }
    delete[] byteArray;
    std::cout << std::endl;
}

int CompressedResidue::preprocess(std::vector<AtomCoordinate>& atoms) {
    int success = 0;
    this->rawAtoms = atoms;
    this->compressedBackBone.resize(atoms.size());
    this->compressedSideChain.resize(atoms.size());
    AminoAcid aa;
    this->AAS = aa.AminoAcids();

    // Remove duplicated atoms
    removeAlternativePosition(atoms);

    // Discretize
    // Extract backbone
    this->backbone = filterBackbone(atoms);

    // this->title = this->strTitle.c_str();
    this->lenTitle = this->strTitle.size();

    this->nResidue = this->backbone.size() / 3;
    this->nBackbone = this->backbone.size();
    this->nAtom = atoms.size();
    this->idxResidue = atoms[0].residue_index;
    this->idxAtom = atoms[0].atom_index;
    this->chain = atoms[0].chain.c_str()[0];
    this->firstResidue = getOneLetterCode(atoms[0].residue);
    this->lastResidue = getOneLetterCode(atoms[atoms.size() - 1].residue);

    // Anchor atoms
    this->_setAnchor();

    if (atoms[atoms.size() - 1].atom == "OXT") {
        this->hasOXT = 1;
        this->OXT = atoms[atoms.size() - 1];
        this->OXT_coords = atoms[atoms.size() - 1].coordinate;
    } else {
        this->hasOXT = 0;
        this->OXT = AtomCoordinate();
        this->OXT_coords = std::vector<float>(3, 0.0);
    }

    std::vector<float> backboneTorsion = getTorsionFromXYZ(this->backbone, 1);
    // Split backbone into phi, psi, omega
    // Calculate phi, psi, omega
    for (int i = 0; i < backboneTorsion.size(); i += 3) {
        this->psi.push_back(backboneTorsion[i]);
        this->omega.push_back(backboneTorsion[i + 1]);
        this->phi.push_back(backboneTorsion[i + 2]);
    }

    std::vector<float> backboneBondAngles = this->nerf.getBondAngles(backbone);
    // Split bond angles into three parts
    for (int i = 1; i < backboneBondAngles.size(); i++) {
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
    this->omegaDisc = Discretizer(this->omega, pow(2, NUM_BITS_OMEGA) - 1);
    this->omegaDiscretized = this->omegaDisc.discretize(this->omega);
    this->psiDisc = Discretizer(this->psi, pow(2, NUM_BITS_PHI_PSI) - 1);
    this->psiDiscretized = this->psiDisc.discretize(this->psi);
    this->n_ca_c_angleDisc = Discretizer(this->n_ca_c_angle, pow(2, NUM_BITS_BOND) - 1);
    this->n_ca_c_angleDiscretized = this->n_ca_c_angleDisc.discretize(this->n_ca_c_angle);
    this->ca_c_n_angleDisc = Discretizer(this->ca_c_n_angle, pow(2, NUM_BITS_BOND) - 1);
    this->ca_c_n_angleDiscretized = this->ca_c_n_angleDisc.discretize(this->ca_c_n_angle);
    this->c_n_ca_angleDisc = Discretizer(this->c_n_ca_angle, pow(2, NUM_BITS_BOND) - 1);
    this->c_n_ca_angleDiscretized = this->c_n_ca_angleDisc.discretize(this->c_n_ca_angle);

    // Set residue names
    this->residueThreeLetter = getResidueNameVector(atoms);

    // Set Discretizer for side chain
    this->sideChainDiscMap = initializeSideChainDiscMap();
    // Calculate side chain info

    this->sideChainAnglesPerResidue = calculateSideChainTorsionAnglesPerResidue(atoms, this->AAS);

    // Discretize side chain
    //this->_discretizeSideChainTorsionAngles(this->sideChainAnglesPerResidue, this->sideChainAnglesDiscretized);
    FixedAngleDiscretizer sideChainDiscretizer = FixedAngleDiscretizer(pow(2, NUM_BITS_TEMP) - 1);
    for (int i = 0; i < this->sideChainAnglesPerResidue.size(); i++) {
        for (int j = 0; j < this->sideChainAnglesPerResidue[i].size(); j++) {
            unsigned int temp = sideChainDiscretizer.discretize(this->sideChainAnglesPerResidue[i][j]);
            this->sideChainAnglesDiscretized.push_back(temp);
        }
    }
    this->nSideChainTorsion = this->sideChainAnglesDiscretized.size();

    // Get tempFactor
    for (int i = 0; i < atoms.size(); i++) {
        this->tempFactors.push_back(atoms[i].tempFactor);
    }
    // Discretize
    this->tempFactorsDisc = Discretizer(this->tempFactors, pow(2, NUM_BITS_TEMP) - 1);
    this->tempFactorsDiscretized = this->tempFactorsDisc.discretize(this->tempFactors);

    // Get header
    this->header = this->get_header();

    // Identify breaks
    this->backboneBreaks = this->nerf.identifyBreaks(backbone);

    // Mark as processed
    this->isPreprocessed = true;

    return success;
}


std::vector<BackboneChain> CompressedResidue::compress(
    std::vector<AtomCoordinate>& atoms
) {
    std::vector<BackboneChain> output;
    // TODO: convert the atom coordinate vector into a vector of compressed residue
    // CURRENT VERSION - 2022-01-10 15:34:21
    // IGNORE BREAKS
    if (!this->isPreprocessed) {
        this->preprocess(atoms);
    }
    this->prevAtoms = {atoms[0], atoms[1], atoms[2]};

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
    res.residue = convertOneLetterCodeToInt(currResCode);
    res.psi = 0; res.omega = 0; res.phi = 0;
    res.n_ca_c_angle = 0; res.ca_c_n_angle = 0; res.c_n_ca_angle = 0;
    output.push_back(res);
    this->compressedBackBone = output;

    this->isCompressed = true;
    return output;
}

std::vector< std::vector<float> > CompressedResidue::_calculateCoordinates() {

    // Iterate through the compressed backbone vector
    std::vector< std::vector<float> > output;
    int total = this->compressedBackBone.size();

    // NOT COMPLETED
    //
    return output;
}

int _restoreResidueNames(
    std::vector<BackboneChain>& compressedBackbone,
    CompressedFileHeader& header,
    std::vector<std::string>& residueThreeLetter
) {
    int success = 0;
    std::string threeLetterCode;
    residueThreeLetter.clear();
    for (int i = 0; i < compressedBackbone.size(); i++) {
        threeLetterCode = convertIntToThreeLetterCode(compressedBackbone[i].residue);
        residueThreeLetter.push_back(threeLetterCode);
    }
    return success;
}

int CompressedResidue::_restoreDiscretizer(int angleType) {
    int success = 0;
    std::vector<unsigned int> temp;

    for (auto cb : this->compressedBackBone) {
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
int CompressedResidue::_restoreAtomCoordinate(float* coords) {
    int success = 0;
    std::string firstResidue = getThreeLetterCode(this->header.firstResidue);
    // TODO: Fix the chain alphabet according to the real data
    // 2022-03-04 14:30:58
    // IMPORTANT: WARNING: TODO: Atom indexing is not correct right now
    // SIDE CHAIN ATOMS SHOULD GE CONSIDERED WHEN INDEXING
    AtomCoordinate prevN = AtomCoordinate(
        "N", firstResidue, "A", this->header.idxAtom, this->header.idxResidue,
        coords[0], coords[1], coords[2]
    );
    AtomCoordinate prevCA = AtomCoordinate(
        "CA", firstResidue, "A", this->header.idxAtom + 1, this->header.idxResidue,
        coords[3], coords[4], coords[5]
    );
    AtomCoordinate prevC = AtomCoordinate(
        "C", firstResidue, "A", this->header.idxAtom + 2, this->header.idxResidue,
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

int CompressedResidue::_getAnchorNum(int threshold) {
    int nAnchor = 0;
    nAnchor = this->nResidue / threshold;
    return nAnchor;
}

void CompressedResidue::_setAnchor() {
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
    std::vector<int> anchorResidueIndices;
    for (int i = 0; i < this->anchorIndices.size(); i++) {
        anchorResidueIndices.push_back(this->anchorIndices[i] + this->idxResidue);
    }
    this->anchorAtoms = getAtomsWithResidueIndex(this->rawAtoms, anchorResidueIndices);
}

std::vector<float> CompressedResidue::checkTorsionReconstruction() {
    std::vector<float> output;
    // Continuize torsion angles
    this->phi = this->phiDisc.continuize(this->phiDiscretized);
    this->psi = this->psiDisc.continuize(this->psiDiscretized);
    this->omega = this->omegaDisc.continuize(this->omegaDiscretized);
    // Append psi, omega, phi to torsion angles
    for (int i = 0; i < this->phi.size(); i++) {
        output.push_back(this->psi[i]);
        output.push_back(this->omega[i]);
        output.push_back(this->phi[i]);
    }
    return output;
}

int CompressedResidue::decompress(std::vector<AtomCoordinate>& atom) {
    // TODO: convert the compressed residue back into a vector of atom coordinate
    // CURRENT VERSION - 2022-02-24 19:12:51
    std::vector<float> torsion_angles;
    std::vector<float> bond_angles;
    int success;

    // Continuize torsion angles
    this->phi = this->phiDisc.continuize(this->phiDiscretized);
    this->psi = this->psiDisc.continuize(this->psiDiscretized);
    this->omega = this->omegaDisc.continuize(this->omegaDiscretized);

    // Append psi, omega, phi to torsion angles
    for (int i = 0; i < (this->phi.size() - 1); i++) {
        torsion_angles.push_back(this->psi[i]);
        torsion_angles.push_back(this->omega[i]);
        torsion_angles.push_back(this->phi[i]);
    }
    // Continuize bond angles
    this->n_ca_c_angle = this->n_ca_c_angleDisc.continuize(this->n_ca_c_angleDiscretized);
    this->ca_c_n_angle = this->ca_c_n_angleDisc.continuize(this->ca_c_n_angleDiscretized);
    this->c_n_ca_angle = this->c_n_ca_angleDisc.continuize(this->c_n_ca_angleDiscretized);
    // Append n_ca_c_angle, ca_c_n_angle, c_n_ca_angle to bond angles
    for (int i = 0; i < this->n_ca_c_angle.size(); i++) {
        bond_angles.push_back(this->ca_c_n_angle[i]);
        bond_angles.push_back(this->c_n_ca_angle[i]);
        bond_angles.push_back(this->n_ca_c_angle[i]);
    }

    // Get the residue vector
    // TODO: EXTRACT RESIDUE NAMES FROM COMPRESSED BACKBONE
    success = _restoreResidueNames(this->compressedBackBone, this->header, this->residueThreeLetter);

    //TODO: FILL IN THIS FUNCTION
    // nerf.reconstruct();
    std::vector<AtomCoordinate> prevForAnchor;
    std::vector<AtomCoordinate> atomByAnchor;
    for (int i = 0; i < this->nAllAnchor - 1; i++) {
        if (i == 0) {
            prevForAnchor = this->prevAtoms;
        }
        // std::vector<int>   sub(&data[100000],&data[101000]);
        std::vector<BackboneChain> subBackbone(
            &this->compressedBackBone[this->anchorIndices[i]],
            &this->compressedBackBone[this->anchorIndices[i + 1]] + 1
        );
        atomByAnchor = reconstructBackboneAtoms(prevForAnchor, subBackbone, this->header);
        // Subset torsion_angles
        std::vector<float> subTorsionAngles(
            &torsion_angles[this->anchorIndices[i] * 3],
            &torsion_angles[this->anchorIndices[i + 1] * 3]
        );
        success = reconstructBackboneReverse(
            atomByAnchor, this->anchorCoordinates[i], subTorsionAngles, this->nerf
        );
        // Append atomByAnchor to atom
        if (i != this->nAllAnchor - 2) {
            atom.insert(atom.end(), atomByAnchor.begin(), atomByAnchor.end() - 3);
        } else {
            atom.insert(atom.end(), atomByAnchor.begin(), atomByAnchor.end());
        }
        // Update prevForAnchor - last 3 atoms of atomByAnchor
        prevForAnchor = std::vector<AtomCoordinate>(
            &atomByAnchor[atomByAnchor.size() - 3],
            &atomByAnchor[atomByAnchor.size()]
        );
    }

    // Reconstruct sidechain
    std::vector< std::vector<AtomCoordinate> > backBonePerResidue = splitAtomByResidue(atom);
    std::string currResidue = getThreeLetterCode(this->header.firstResidue);
    AminoAcid aa;
    std::map<std::string, AminoAcid> aminoAcidMap = aa.AminoAcids();
    this->AAS =aminoAcidMap;
    std::vector<AtomCoordinate> fullResidue;

    success = this->_continuizeSideChainTorsionAngles(
        this->sideChainAnglesDiscretized, this->sideChainAnglesPerResidue
    );

    for (int i = 0; i < backBonePerResidue.size(); i++) {
        if (i != 0){
            currResidue = backBonePerResidue[i][0].residue;
        }
        fullResidue = nerf.reconstructAminoAcid(
            backBonePerResidue[i], this->sideChainAnglesPerResidue[i], aminoAcidMap[currResidue]
        );
        if (this->useAltAtomOrder) {
            _reorderAtoms(fullResidue, aminoAcidMap[currResidue]);
        }
        backBonePerResidue[i] = fullResidue;
    }
    // Flatten backBonePerResidue
    atom.clear();
    for (int i = 0; i < backBonePerResidue.size(); i++) {
        for (int j = 0; j < backBonePerResidue[i].size(); j++) {
            atom.push_back(backBonePerResidue[i][j]);
        }
    }
    // Reindex atom index of atom
    if (this->hasOXT) {
        atom.push_back(this->OXT);
    }
    setAtomIndexSequentially(atom, this->header.idxAtom);

    // Set tempFactor
    std::vector<float> tempFactors = this->tempFactorsDisc.continuize(this->tempFactorsDiscretized);
    for (int i = 0; i < atom.size(); i++) {
        atom[i].tempFactor = tempFactors[i];
    }

    this->rawAtoms = atom;

    return success;
}

int CompressedResidue::read(std::string filename) {
    int success;
    // Open file in reading binary mode
    std::ifstream file(filename, std::ios::binary);
    // Check if file is open
    if (!file.is_open()) {
        std::cout << "Error: Could not open file " << filename << std::endl;
        return -1;
    }
    // Check the file starts with magic number
    char mNum[MAGICNUMBER_LENGTH];
    file.read(mNum, MAGICNUMBER_LENGTH);
    // compare mNum and this->magicNumber
    for (int i = 0; i < MAGICNUMBER_LENGTH; i++) {
        if (mNum[i] != MAGICNUMBER[i]) {
            std::cout << "Error: File " << filename << " is not a valid compressed residue file" << std::endl;
            return -1;
        }
    }
    // Read the header
    file.read(reinterpret_cast<char*>(&this->header), sizeof(this->header));
    this->read_header(this->header);
    // Read anchorIndices
    this->anchorIndices.resize(this->nAllAnchor);
    // Read int vector
    file.read(reinterpret_cast<char*>(&this->anchorIndices[0]), sizeof(int) * this->nAllAnchor);
    // Read the title
    char title[this->header.lenTitle];
    file.read(title, sizeof(char) * this->header.lenTitle);
    this->strTitle = std::string(title, (int)this->header.lenTitle);
    // Read the prev atoms
    // In the file, only xyz coordinates are stored
    // So, we need to reconstruct the atomcoordinate from the xyz coordinates & the information from the header
    float prevAtomCoords[9];
    file.read(reinterpret_cast<char*>(prevAtomCoords), sizeof(prevAtomCoords));

    // ANCHOR ATOMS
    if (this->header.nAnchor > 2) {
        float innerAnchorCoords[9 * (this->header.nAnchor - 2)];
        file.read(reinterpret_cast<char*>(innerAnchorCoords), sizeof(innerAnchorCoords));
        for (int i = 0; i < (this->header.nAnchor - 2); i++) {
            std::vector< std::vector<float> > innerAnchorCoord;
            for (int j = 0; j < 3; j++) {
                innerAnchorCoord.push_back({
                    innerAnchorCoords[i * 9 + j * 3],
                    innerAnchorCoords[i * 9 + j * 3 + 1],
                    innerAnchorCoords[i * 9 + j * 3 + 2]
                });
            }
            this->anchorCoordinates.push_back(innerAnchorCoord);
        }
    }

    // LAST ATOMS
    float lastAtomCoords[9];
    file.read(reinterpret_cast<char*>(lastAtomCoords), sizeof(lastAtomCoords));
    for (int i = 0; i < 3; i++) {
        this->lastAtomCoordinates.push_back({ lastAtomCoords[i*3], lastAtomCoords[i*3 + 1], lastAtomCoords[i*3+ 2] });
    }
    this->anchorCoordinates.push_back(this->lastAtomCoordinates);

    file.read(&this->hasOXT, sizeof(char));
    float oxtCoords[3];
    file.read(reinterpret_cast<char*>(oxtCoords), sizeof(oxtCoords));
    this->OXT_coords = { oxtCoords[0], oxtCoords[1], oxtCoords[2] };
    this->OXT = AtomCoordinate(
        "OXT", getThreeLetterCode(this->header.lastResidue), std::string(1, this->header.chain),
        (int)this->header.nAtom, (int)this->header.nResidue, this->OXT_coords
    );

    // Read sidechain discretizers
    // file.read(reinterpret_cast<char*>(&this->sideChainDisc), sizeof(SideChainDiscretizers));time

    // Read the array of backbone bytes
    // Read 8 bytes at a time
    // Converted data will be saved at this->compressedBackBone
    // if this->compressedBackBone is not empty, then it will be overwritten
    // check if compressedBackBone is empty
    if (this->compressedBackBone.size() == 0) {
        this->compressedBackBone.resize(this->header.nResidue);
    } else {
        this->compressedBackBone.clear();
        this->compressedBackBone.resize(this->header.nResidue);
    }
    char* buffer = new char[8];
    for (int i = 0; i < this->header.nResidue; i++) {
        file.read(buffer, 8);
        // Convert buffer to compressedBackBone using convertBytesToBackboneChain
        this->compressedBackBone[i] = convertBytesToBackboneChain(buffer);
    }
    delete[] buffer;

    // Read sidechain
    int encodedSideChainSize = this->header.nSideChainTorsion;
    if (encodedSideChainSize % 2 == 1) {
        encodedSideChainSize++;
    }
    encodedSideChainSize /= 2;
    unsigned char* encodedSideChain = new unsigned char[this->header.nSideChainTorsion];
    // file.read(reinterpret_cast<char*>(encodedSideChain), sizeof(encodedSideChain));
    // Decode sidechain
    // success = decodeSideChainTorsionVector(
    //     encodedSideChain, this->header.nSideChainTorsion,
    //     this->sideChainAnglesDiscretized
    // );
    file.read(reinterpret_cast<char*>(encodedSideChain), this->header.nSideChainTorsion);
    // Read char array to unsigned int vector sideChainAnglesDiscretized
    unsigned int temp;
    for (int i = 0; i < this->header.nSideChainTorsion; i++) {
        // Convert char to unsigned int
        temp = encodedSideChain[i];
        this->sideChainAnglesDiscretized.push_back(temp);
    }
    delete[] encodedSideChain;

    // Read temperature factor
    // Read discretizer for temperature factor
    file.read(reinterpret_cast<char*>(&this->tempFactorsDisc.min), sizeof(float));
    file.read(reinterpret_cast<char*>(&this->tempFactorsDisc.cont_f), sizeof(float));

    unsigned char* encodedTempFactors = new unsigned char[this->header.nAtom];
    file.read(reinterpret_cast<char*>(encodedTempFactors), this->header.nAtom);

    unsigned int tempFactor;
    for (int i = 0; i < this->header.nAtom; i++) {
        tempFactor = (unsigned int)encodedTempFactors[i];
        this->tempFactorsDiscretized.push_back(tempFactor);
    }
    delete[] encodedTempFactors;

    success = _restoreAtomCoordinate(prevAtomCoords);
    if (success != 0) {
        std::cout << "Error: Could not restore prevAtoms" << std::endl;
        return -1;
    }
    for (int i = 0; i < 6; i++) {
        success = _restoreDiscretizer(i);
    }

    // Close file
    file.close();
    return success;
}

int CompressedResidue::write(std::string filename) {
    int flag = 0;
    std::ofstream outfile;
    outfile.open(filename, std::ios::out | std::ios::binary);
    // Open in binary & writing mode
    if (!outfile.is_open()) {
        std::cout << "Error opening file: " << filename << std::endl;
        flag = -1;
        return flag;
    }
    if (outfile.good()) {
        // Write magic number
        outfile.write(MAGICNUMBER, MAGICNUMBER_LENGTH);
        // Write header
        outfile.write((char*)&this->header, sizeof(CompressedFileHeader));
        // Write anchorIndices
        for (int i = 0; i < this->anchorIndices.size(); i++) {
            outfile.write((char*)&this->anchorIndices[i], sizeof(int));
        }
        // Write title
        char* title = new char[this->strTitle.length() + 1];
        strcpy(title, this->strTitle.c_str());
        outfile.write((char*)title, sizeof(char) * this->header.lenTitle);
        delete[] title;

        // 2022-08-08 19:15:30 - Changed to write all anchor atoms
        // TODO: NEED TO BE CHECKED
        for (auto anchors : this->anchorAtoms) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    outfile.write((char*)&anchors[i].coordinate[j], sizeof(float));
                }
            }
        }

        outfile.write(&this->hasOXT, sizeof(char));
        for (int i = 0; i < 3; i++) {
            outfile.write((char*)&this->OXT_coords[i], sizeof(float));
        }

        // END OF BACKBONE METADATA
        // // Write sidechain discretizers
        // outfile.write((char*)&this->sideChainDisc, sizeof(SideChainDiscretizers));

        // Write the compressed backbone
        char* buffer = new char[8];
        for (int i = 0; i < this->compressedBackBone.size(); i++) {
            flag = convertBackboneChainToBytes(this->compressedBackBone[i], buffer);
            outfile.write(buffer, 8);
        }
        delete[] buffer;
        // 2022-06-07 18:44:09 - Check memory leak
        // Write side chain torsion angles
        int encodedSideChainSize = this->nSideChainTorsion;
        if (encodedSideChainSize % 2 == 1) {
            encodedSideChainSize++;
        }
        encodedSideChainSize /= 2;
        // char* charSideChainTorsion = encodeSideChainTorsionVector(this->sideChainAnglesDiscretized);
        // outfile.write(charSideChainTorsion, encodedSideChainSize);
        // TODO: TESTING
        // Convert unsigned int to char array
        unsigned char* charSideChainTorsion = new unsigned char[this->nSideChainTorsion];
        // Get array of unsigned int from sideChainAnglesDiscretized and convert to char array
        for (int i = 0; i < this->nSideChainTorsion; i++) {
            // convert unsigned int to char
            charSideChainTorsion[i] = (unsigned char)this->sideChainAnglesDiscretized[i];
        }
        outfile.write((char*)charSideChainTorsion, this->sideChainAnglesDiscretized.size());
        delete[] charSideChainTorsion;

        // Write temperature factors
        // Disc
        outfile.write((char*)&this->tempFactorsDisc.min, sizeof(float));
        outfile.write((char*)&this->tempFactorsDisc.cont_f, sizeof(float));

        unsigned char* charTempFactors = new unsigned char[this->header.nAtom];
        for (int i = 0; i < this->header.nAtom; i++) {
            charTempFactors[i] = (unsigned char)this->tempFactorsDiscretized[i];
        }
        outfile.write((char*)charTempFactors, this->tempFactorsDiscretized.size());
        delete[] charTempFactors;

        // Close file
        outfile.close();
    } else {
        std::cout << "Error writing file: " << filename << std::endl;
        flag = -1;
    }
    return flag;
}

CompressedFileHeader CompressedResidue::get_header() {
    CompressedFileHeader header;
    // counts
    header.nResidue = this->nResidue;
    header.nAtom = this->nAtom;
    header.idxResidue = this->idxResidue;
    header.idxAtom = this->idxAtom;
    header.nAnchor = this->nAllAnchor;
    header.nSideChainTorsion = this->nSideChainTorsion;
    header.firstResidue = this->firstResidue;
    header.lastResidue = this->lastResidue;
    header.lenTitle = this->lenTitle;
    //
    header.chain = this->chain;
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

int CompressedResidue::read_header(CompressedFileHeader& header) {
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

void CompressedResidue::print(int length) {
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
        std::cout << "x: " << this->prevAtoms[i].coordinate[0] << std::endl;
        std::cout << "y: " << this->prevAtoms[i].coordinate[1] << std::endl;
        std::cout << "z: " << this->prevAtoms[i].coordinate[2] << std::endl;
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

void CompressedResidue::printSideChainTorsion(std::string filename) {
    std::ofstream outfile;
    outfile.open(filename);
    //
    outfile << "ResidueInd,Residue,Type,Key,RawVal,DiscVal,DiscMin,DiscContF,ReconVal,Diff\n";
    int movingIndex = 0;
    std::vector< std::vector<AtomCoordinate> > atomByResidue = splitAtomByResidue(this->rawAtoms);
    std::map<std::string, float> currBondAngle;
    std::map<std::string, float> currBondLength;
    std::map<std::string, float> currTorsionAngle;

    for (int i = 0; i < this->nResidue; i++) {
        currBondAngle = calculateBondAngles(atomByResidue[i], this->AAS[this->residueThreeLetter[i]]);
        currBondLength = calculateBondLengths(atomByResidue[i], this->AAS[this->residueThreeLetter[i]]);
        for (auto bl: currBondLength) {
            outfile << i << "," << this->residueThreeLetter[i] << ",BondLength,";
            outfile << bl.first << "," << bl.second << ",NA,NA,NA,";
            outfile << this->AAS[this->residueThreeLetter[i]].bondLengths[bl.first] << ",";
            outfile << bl.second - this->AAS[this->residueThreeLetter[i]].bondLengths[bl.first] << "\n";
        }
        for (auto ba : currBondAngle) {
            outfile << i << "," << this->residueThreeLetter[i] << ",BondAngle,";
            outfile << ba.first << "," << ba.second << ",NA,NA,NA," << this->AAS[this->residueThreeLetter[i]].bondAngles[ba.first] << ",";
            outfile << ba.second - this->AAS[this->residueThreeLetter[i]].bondAngles[ba.first] << "\n";
        }
        for (int j = 0; j < this->sideChainAnglesPerResidue[i].size(); j++) {
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

void _reorderAtoms(std::vector<AtomCoordinate>& atoms, AminoAcid& aa) {
    std::vector<AtomCoordinate> newAtoms = atoms;
    for (int i = 0; i < atoms.size(); i++) {
        if (atoms[i].atom == aa.altAtoms[i]) {
            continue;
        } else {
            for (int j = 0; j < aa.altAtoms.size(); j++) {
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
    for (int i = 0; i < newSize; i++) {
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
    for (int i = 0; i < vector.size(); i++) {
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
    for (int i = 0; i < aaNames.size(); i++) {
        std::vector<Discretizer> disc;
        numTorsion = getSideChainTorsionNum(aaNames[i]);
        for (int j = 0; j < numTorsion; j++) {
            Discretizer currDisc = Discretizer();
            disc.push_back(currDisc);
        }
        discMap[aaNames[i]] = disc;
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

