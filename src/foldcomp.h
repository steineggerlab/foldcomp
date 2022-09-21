/**
 * File: foldcomp.h
 * Project: foldcomp
 * Created: 2021-02-02 14:04:40
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     This file contains main data structures for torsion angle compression and
 *     functions for handling them.
 * ---
 * Last Modified: 2022-09-21 19:56:58
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
#pragma once
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <bitset>
#include <string.h>
#include "amino_acid.h"
#include "atom_coordinate.h"
#include "discretizer.h"
#include "sidechain.h"
#include "torsion_angle.h"
#include "nerf.h"
#include "utility.h"

// TAR format handling
#include "microtar/microtar.h"

// CONSTANTS
#define NUM_TYPE_OF_ANGLES 6
#define MAGICNUMBER_LENGTH 4
#define MAGICNUMBER "FCMP"
// 2022-04-14 15:37:51
#define NUM_BITS_PHI_PSI 12
#define NUM_BITS_OMEGA 11
#define NUM_BITS_BOND 8
#define NUM_BITS_RESIDUE 5
#define NUM_BITS_TEMP 8
#define NUM_BITS_SIDECHAIN 4

#define N_TO_CA_DIST 1.4581
#define CA_TO_C_DIST 1.5281
#define C_TO_N_DIST 1.3311
#define PRO_N_TO_CA_DIST 1.353

enum ValidityError {
    SUCCESS = 0,
    E_WRONG_MAGIC_NUMBER,
    E_BACKBONE_COUNT_MISMATCH,
    E_SIDECHAIN_COUNT_MISMATCH,
    E_TEMP_FACTOR_COUNT_MISMATCH,
    E_EMPTY_BACKBONE_ANGLE,
    E_EMPTY_SIDECHAIN_ANGLE,
    E_EMPTY_TEMP_FACTOR
};

// NOTE: ORDER OF BOND ANGLE: CA_C_N, C_N_CA, N_CA_C
// NOTE: THE ORDER OF TORSION ANGLE IS PSI->OMEGA->PHI
struct BackboneChain {
    /* data */
    // TOTAL BITS: 64 bits = 8 bytes
    unsigned int residue: NUM_BITS_RESIDUE; // 5 bits
    unsigned int omega : NUM_BITS_OMEGA;   // 11 bits
    unsigned int psi : NUM_BITS_PHI_PSI;    // 12 bits
    unsigned int phi : NUM_BITS_PHI_PSI;    // 12 bits
    unsigned int ca_c_n_angle : NUM_BITS_BOND;  // 8 bits
    unsigned int c_n_ca_angle : NUM_BITS_BOND;  // 8 bits
    unsigned int n_ca_c_angle: NUM_BITS_BOND;  // 8 bits
};

// IDEA: Split backbone chain header & general header??
// First residue string should be saved in the header
// TODO: SAVE FIRST RESIDUE
struct BackboneChainHeader {
    unsigned int nResidue: 16; // 16 bits
    unsigned int nAtom: 16;    // 16 bits
    unsigned int idxResidue: 16; // 16 bits
    unsigned int idxAtom: 16;    // 16 bits
    // char firstResidue; // 1 byte = 8 bits --> WILL BE APPLIED AFTER CHANGING BACKBONE CHAIN
    float prevAtoms[9];
};

struct DecompressedBackboneChain {
    /* data */
    char residue;
    float n_ca_c_angle;
    float ca_c_n_angle;
    float c_n_ca_angle;
    float phi;
    float psi;
    float omega;
};

// 2022-06-16 21:35:20 - REMOVED UNNECESSARY STRUCTURES

// TODO: Implement offset array
// An array of offsets?
struct FileOffset {
    unsigned int firstBackboneStart: 32;
    unsigned int sideChainStart: 32;
};


// TODO: Split fixed size header & variable size header
// An array of chain names


struct CompressedFileHeader {
    /* data */
    // TOTAL:
    // Backbone
    unsigned int nResidue: 16; // 16 bits
    unsigned int nAtom: 16;    // 16 bits
    unsigned int idxResidue: 16; // 16 bits
    unsigned int idxAtom: 16;    // 16 bits
    // Anchor points - 2022-08-08 18:28:51
    unsigned int nAnchor: 8; // 8 bits
    char chain;
    // Sidechain
    unsigned int nSideChainTorsion: 32;
    char firstResidue;
    char lastResidue;
    unsigned int lenTitle: 32;
    // Discretizer for backbone chain
    float mins[6];
    float cont_fs[6];
    // IDEA: Offset?
};

struct CompressedFileOffset {
    unsigned int firstBackboneStart: 32;
    unsigned int sideChainStart: 32;
};

struct SideChainDiscretizers {
    float ala_min[2];
    float ala_cont_fs[2];
    float arg_min[8];
    float arg_cont_fs[8];
    float asn_min[5];
    float asn_cont_fs[5];
    float asp_min[5];
    float asp_cont_fs[5];
    float cys_min[3];
    float cys_cont_fs[3];
    float gln_min[6];
    float gln_cont_fs[6];
    float glu_min[6];
    float glu_cont_fs[6];
    float gly_min[1];
    float gly_cont_fs[1];
    float his_min[7];
    float his_cont_fs[7];
    float ile_min[5];
    float ile_cont_fs[5];
    float leu_min[5];
    float leu_cont_fs[5];
    float lys_min[6];
    float lys_cont_fs[6];
    float met_min[5];
    float met_cont_fs[5];
    float phe_min[6];
    float phe_cont_fs[6];
    float pro_min[4];
    float pro_cont_fs[4];
    float ser_min[3];
    float ser_cont_fs[3];
    float thr_min[4];
    float thr_cont_fs[4];
    float trp_min[11];
    float trp_cont_fs[11];
    float tyr_min[9];
    float tyr_cont_fs[9];
    float val_min[4];
    float val_cont_fs[4];
};

struct SidechainAngles {
    unsigned int torsion1: 4;
    unsigned int torsion2: 4;
};

// Conversion
uint32_t convertCompressedResidueToFirst4Bytes(BackboneChain& res);
uint32_t convertCompressedResidueToSecond4Bytes(BackboneChain& res);
int convertBackboneChainToBytes(BackboneChain& res, char* output);
BackboneChain convertBytesToBackboneChain(char* bytes);
BackboneChain newBackboneChain(
    char residue, unsigned int phi, unsigned int psi, unsigned int omega,
    unsigned int n_ca_c_angle, unsigned int ca_c_n_angle, unsigned int c_n_ca_angle
);
BackboneChain newBackboneChain(
    unsigned int bResidue, unsigned int phi, unsigned int psi, unsigned int omega,
    unsigned int n_ca_c_angle, unsigned int ca_c_n_angle, unsigned int c_n_ca_angle
);

// TODO: Change header
DecompressedBackboneChain decompressBackboneChain(
    BackboneChain& bb, CompressedFileHeader& header
);
std::vector<DecompressedBackboneChain> decompressBackboneChain(
    std::vector<BackboneChain>& bbv, CompressedFileHeader& header
);

float _continuize(unsigned int input, float min, float cont_f);

// Reconstruct
std::vector<AtomCoordinate> reconstructBackboneAtoms(
    std::vector<AtomCoordinate>& prevAtoms,
    std::vector<BackboneChain>& backbone,
    CompressedFileHeader& header
);

std::vector<AtomCoordinate> reconstructSidechainAtoms(
    std::vector<AtomCoordinate>& backBoneAtoms,
    std::vector<AtomCoordinate>& wholeAtoms,
    SideChainDiscretizers& sideChainDisc
);

int discretizeSideChainTorsionAngles(
    std::vector< std::vector<float> >& torsionPerResidue,
    std::vector<std::string>& residueNames,
    std::map<std::string, AminoAcid>& AAS,
    SideChainDiscretizers& scDiscretizers,
    std::map<std::string, std::vector<Discretizer> >& scDiscretizersMap,
    std::vector<unsigned int>& output
);


int continuizeSideChainTorsionAngles(
    std::vector<unsigned int>& torsionDiscretized,
    std::vector<std::string>& residueNames,
    std::map<std::string, AminoAcid>& AAS,
    SideChainDiscretizers& scDiscretizers,
    std::map<std::string, std::vector<Discretizer> >& scDiscretizersMap,
    std::vector< std::vector<float> >& output
);

float* getContFFromSideChainDiscretizers(
    std::string& residue, SideChainDiscretizers& scDiscretizers
);
float* getMinPointerFromSideChainDiscretizers(
    std::string& residue, SideChainDiscretizers& scDiscretizers
);
std::map<std::string, std::vector<Discretizer> > initializeSideChainDiscMap();
unsigned char* encodeDiscretizedTempFactors(std::vector<unsigned int> vector);
int decodeDiscretizedTempFactors(unsigned char* input, int size, std::vector<unsigned int>& vector);
char* encodeSideChainTorsionVector(std::vector<unsigned int> vector);
int decodeSideChainTorsionVector(char* input, int nTorsion, std::vector<unsigned int>& vector);
int getSideChainTorsionNum(std::string residue);
int fillSideChainDiscretizerMap(
    SideChainDiscretizers& scDiscretizers,
    std::map<std::string, std::vector<Discretizer> >& scDiscretizersMap
);

void _reorderAtoms(std::vector<AtomCoordinate>& atoms, AminoAcid& aa);

// Print
void printCompressedResidue(BackboneChain& res);
void printValidityError(ValidityError err, std::string& filename);

struct FloatArrayWithDisc {
    unsigned short size;
    float min;
    float cont_f;
    float* array;
};

class Foldcomp {
private:
    /* data */
    std::string magicNumber;
    /* private methods */
    std::vector< std::vector<float> > _calculateCoordinates();
    int _restoreDiscretizer(int angleType);
    int _restoreAtomCoordinate(float* coords);
    int _preprocessBackbone();
    int _preprocessSideChain();

    // Anchor
    int _getAnchorNum(int threshold);
    void _setAnchor();
    std::vector<AtomCoordinate> _getAnchorAtoms(bool includeStartAndEnd = true);

    int _discretizeSideChainTorsionAngles(
        std::vector< std::vector<float> >& input, std::vector<unsigned int>& output
    );

    int _continuizeSideChainTorsionAngles(
        std::vector<unsigned int>& input, std::vector< std::vector<float> >& output
    );
    // Check validity
    ValidityError _checkValidity();
public:
    Foldcomp(/* args */){};
    ~Foldcomp(){};
    bool isPreprocessed = false;
    bool isCompressed = false;
    bool backwardReconstruction = true;
    bool useAltAtomOrder = false;
    bool checkValidity = true;
    // Number of atoms & residues
    int nResidue = 0;
    int nAtom = 0;
    int nBackbone = 0;
    int nSideChainTorsion = 0;
    int nInnerAnchor = 0;
    int nAllAnchor = 0;
    int anchorThreshold = 200;
    // Indices for residue & atom
    int idxResidue = 0;
    int idxAtom = 0;
    char chain;
    char firstResidue;
    char lastResidue;
    char hasOXT = 1;

    // Metadata
    std::string strTitle;
    // const char* title;
    int lenTitle;
    std::map<std::string, std::string> strMetadata;
    std::map<std::string, std::vector<float> > floatMetadata;

    // Header
    CompressedFileHeader header;

    // Vectors for atoms
    std::vector<AtomCoordinate> rawAtoms;
    std::vector<AtomCoordinate> prevAtoms; // 3 atoms
    std::vector<AtomCoordinate> lastAtoms; // 3 atoms
    std::vector< std::vector<float> > lastAtomCoordinates; // 3 atoms
    std::vector<AtomCoordinate> backbone;
    // 2022-08-05 20:47:57 - Anchors for reducing RMSD
    std::vector< std::vector<AtomCoordinate> > anchorAtoms;
    std::vector< std::vector< std::vector<float> > > anchorCoordinates;
    std::vector<int> anchorIndices;

    std::vector<BackboneChain> compressedBackBone;
    std::vector<unsigned int> compressedSideChain;
    std::vector<int> backboneBreaks;
    std::vector<char> residues;
    std::vector<std::string> residueThreeLetter;
    AtomCoordinate OXT;
    std::vector<float> OXT_coords;

    // Angles
    std::vector<float> psi;
    std::vector<float> omega;
    std::vector<float> phi;
    std::vector<float> n_ca_c_angle;
    std::vector<float> ca_c_n_angle;
    std::vector<float> c_n_ca_angle;
    Discretizer psiDisc;
    std::vector<unsigned int> psiDiscretized;
    Discretizer omegaDisc;
    std::vector<unsigned int> omegaDiscretized;
    Discretizer phiDisc;
    std::vector<unsigned int> phiDiscretized;
    Discretizer n_ca_c_angleDisc;
    std::vector<unsigned int> n_ca_c_angleDiscretized;
    Discretizer ca_c_n_angleDisc;
    std::vector<unsigned int> ca_c_n_angleDiscretized;
    Discretizer c_n_ca_angleDisc;
    std::vector<unsigned int> c_n_ca_angleDiscretized;
    std::map<std::string, AminoAcid> AAS;
    Nerf nerf;
    // Sidechain angles
    std::vector<float> sideChainAngles;
    std::vector< std::vector<float> > sideChainAnglesPerResidue;
    std::vector<unsigned int> sideChainAnglesDiscretized;
    SideChainDiscretizers sideChainDisc;
    std::map<std::string, std::vector<Discretizer> > sideChainDiscMap;
    int encodedSideChainSize;
    // Temperature factors / pLDDT
    std::vector<float> tempFactors;
    std::vector<unsigned int> tempFactorsDiscretized;
    Discretizer tempFactorsDisc;

     // methods
    int preprocess(std::vector<AtomCoordinate>& atoms);
    std::vector<BackboneChain> compress(std::vector<AtomCoordinate>& atoms);
    int decompress(std::vector<AtomCoordinate>& atoms);
    int read(std::istream & filename);
    int write(std::string filename);
    // Read & write for tar files
    // int readTar(mtar_t& tar, std::string filename, size_t size);
    int writeTar(mtar_t& tar, std::string filename, size_t size);

    int reconstruct(std::vector<AtomCoordinate>& atoms, int mode);
    CompressedFileHeader get_header();
    int read_header(CompressedFileHeader& header);
    size_t getSize();
    // methods for getting plddt (tempFactors) or amino acid sequence
    int continuizeTempFactors();
    int writeFASTALike(std::string filename, std::vector<std::string>& data);
    int writeFASTALikeTar(mtar_t& tar, std::string filename, std::vector<std::string>& data);
    int extract(std::vector<std::string>& data, int type);

    // temporary method for testing
    std::vector<float> checkTorsionReconstruction();
    void print(int length = 5);
    void printSideChainTorsion(std::string filename);
};
