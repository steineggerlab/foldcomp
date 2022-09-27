/**
 * File: nerf.h
 * Project: foldcomp
 * Created: 2021-01-11 16:42:32
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     This file contains NeRF (Natural Extension of Reference Frame) algorithm
 *     implementation. NeRF is a method to calculate upcoming atoms' coordinates
 *     with three preceding atoms and bond information.
 *     Reference: https://benjamin.computer/posts/2018-03-16-mres-part2.html
 * ---
 * Last Modified: 2022-07-20 01:55:43
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */

// COMMENT FROM REFERENCE
// # TODO - PROLINE has different lengths which we should take into account
// # TODO - A_TO_C angle differs by + / -5 degrees

#pragma once

#include <fstream>
#include <map>
#include <string>
#include <vector>

class AminoAcid;
class AtomCoordinate;

class Nerf {
private:
    // private members
public:
    Nerf(/* args */);
    ~Nerf();
    /* data */
    std::map<std::string, float> bond_lengths {
        // TODO: if residue is proline, apply different bond length
        // 2021-01-19 15:20:42 NOTE: TEST FORE SERINE
        {"N_TO_CA", 1.4581}, {"PRO_N_TO_CA", 1.353}, // proline has different lengths
        {"CA_TO_C", 1.5281}, {"C_TO_N", 1.3311},
        {"C_TO_O", 1.23}
    };
    // bond angles are in radian
    std::map<std::string, float> bond_angles {
        {"N_TO_CA", 121.3822}, {"CA_TO_C", 111.2812}, {"C_TO_N", 116.6429},
        {"C_TO_O", 120.5117}
    };
    //NOTE: bond_order??
    // PDB file is ordered from N terminal to C terminal
    std::vector<std::string> bond_order { "C_TO_N", "N_TO_CA", "CA_TO_C" };

    std::vector<float> place_atom(
        std::vector< std::vector<float> > prev_atoms, // Need three previous atoms
        float bond_length, float bond_angle,
        float torsion_angle
    );

    // 2021-01-22 13:41:03
    // NOTE: Testing 3 cases for reconstruction
    //  1) Use constants for bond angles and bond lengths
    //  2) Use peptide-specific bond angles and bond lengths
    //  3) Use dynamic (atom specific) bond angles and bond lengths

    std::vector<AtomCoordinate> reconstructWithConstants(
        std::vector<AtomCoordinate> original_atoms,
        std::vector<float> torsion_angles);

    std::vector<AtomCoordinate> reconstructWithAASpecificAngles(
        std::vector<AtomCoordinate> original_atoms,
        std::vector<float> torsion_angles,
        std::map<std::string, float> aa_bond_lengths,
        std::map< std::string, std::map<std::string, float> > aa_bond_angles
    );

    std::vector<AtomCoordinate> reconstructWithAAMaps(
        std::vector<AtomCoordinate> original_atoms,
        std::map < std::string, std::vector<std::string> > prev_atom_map,
        std::map<std::string, float> aa_torsion_angles,
        std::map<std::string, float> aa_bond_lengths,
        std::map<std::string, float> aa_bond_angles
    );

    std::vector<AtomCoordinate> reconstructWithDynamicAngles(
        std::vector<AtomCoordinate> original_atoms,
        std::vector<float> torsion_angles,
        std::vector<float> atom_bond_lengths,
        std::vector<float> atom_bond_angles
    );

    std::vector<AtomCoordinate> reconstructWithDynamicAngles(
        std::vector<AtomCoordinate> original_atoms,
        std::vector<float> torsion_angles,
        std::vector<float> atom_bond_angles
    );

    std::vector<AtomCoordinate> reconstructWithBreaks(
        std::vector<AtomCoordinate> original_atoms,
        std::vector<float> torsion_angles,
        std::vector<float> atom_bond_angles,
        std::vector<int> break_indices
    );

    std::vector<AtomCoordinate> reconstrutWithRelativePositions(
        std::vector<AtomCoordinate> original_atoms,
        std::vector<AtomCoordinate> relative_position
    );

    std::vector<AtomCoordinate> reconstructWithReversed(
        std::vector<AtomCoordinate> original_atoms,
        std::vector<float> torsion_angles,
        std::vector<float> atom_bond_angles
    );

    std::vector<AtomCoordinate> reconstructAminoAcid(
        std::vector<AtomCoordinate> original_atoms,
        std::vector<float> torsion_angles,
        AminoAcid& aa
    );

    void writeInfoForChecking(
        std::vector<AtomCoordinate> coord_list,
        std::string output_path, std::string sep = ","
    );
    void writeCoordinatesBinary(
        std::vector<AtomCoordinate> coord_list, std::string output_path
    );
    std::vector<float> getBondAngles(std::vector<AtomCoordinate> original_atoms);
    std::vector<float> getBondLengths(std::vector<AtomCoordinate> original_atoms);

    // 2021-07-12 17:18:16
    // Method to identify breaks from AtomCoordinate vector
    std::vector<int> identifyBreaks(
        std::vector<AtomCoordinate> original_atoms, float cutoff = 2
    );


};
