/**
 * File: nerf.cpp
 * Project: foldcomp
 * Created: 2021-01-11 16:42:27
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     This file contains NeRF (Natural Extension of Reference Frame) algorithm
 *     implementation. NeRF is a method to calculate upcoming atoms' coordinates
 *     with three preceding atoms and bond information.
 *     Reference: https://benjamin.computer/posts/2018-03-16-mres-part2.html
 * ---
 * Last Modified: 2022-09-13 15:15:02
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */

#include "nerf.h"


Nerf::Nerf(/* args */) {}

Nerf::~Nerf() {}

/**
 * @brief Calculate the position of next atom with torsion angle.
 *
 * @details
 * This code is a part of FoldU project. this function calculates the
 * position of the next atom with NeRF algorithm.
 * Angles are in radians, lengths in angstroms
 * @param prev_atoms XYZ coordinates of 3 previous atoms
 * @param bond_length a float, the length of current bond
 * @param bond_angle a float, the angle of current bond
 * @param torsion_angle a float, the torsion angle of current bond.
 * @return std::vector<float>
 */
std::vector<float> Nerf::place_atom(
    std::vector< std::vector<float> > prev_atoms,
    float bond_length, float bond_angle,
    float torsion_angle
){
    // 00. Get 3 atom coordinates
    std::vector<float> atm_a = prev_atoms[0];
    std::vector<float> atm_b = prev_atoms[1];
    std::vector<float> atm_c = prev_atoms[2];

    // 01. Obtain vectors from coordinates

    std::vector<float> ab {
        (atm_b[0] - atm_a[0]), (atm_b[1] - atm_a[1]), (atm_b[2] - atm_a[2])
    };
    std::vector<float> bc {
        (atm_c[0] - atm_b[0]), (atm_c[1] - atm_b[1]), (atm_c[2] - atm_b[2])
    };
    float bc_norm = norm(bc);

    // n2 - unit vector of direction same with d2
    std::vector<float> bcn {
        (bc[0] / bc_norm), (bc[1] / bc_norm), (bc[2] / bc_norm)
    };
    // Current atom
    bond_angle = bond_angle * 3.14159265 / 180.0;
    torsion_angle = torsion_angle * 3.14159265 / 180.0;

    std::vector<float> curr_atm {
        (-1 * bond_length * cosf(bond_angle)), // x
        (bond_length * cosf(torsion_angle) * sinf(bond_angle)), // y
        (bond_length * sinf(torsion_angle) * sinf(bond_angle)) // z
    };

    // 02. Calculate cross product
    // TODO: check if we can reduce the count of norm calculation
    std::vector<float> n = crossProduct(ab, bcn);
    float n_norm = norm(n);
    n[0] = n[0] / n_norm;
    n[1] = n[1] / n_norm;
    n[2] = n[2] / n_norm;
    std::vector<float> nbc = crossProduct(n, bcn);
    std::vector< std::vector<float> > m {
        {bcn[0], nbc[0], n[0]},
        {bcn[1], nbc[1], n[1]},
        {bcn[2], nbc[2], n[2]}
    };
    //
    std::vector<float> atm_d {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            atm_d[i] += (m[i][j] * curr_atm[j]);
        }
    }

    // 0X. Add the coordinates of atom c to atom d --> absolute coordinate
    atm_d[0] += atm_c[0];
    atm_d[1] += atm_c[1];
    atm_d[2] += atm_c[2];

    return atm_d;
}

/**
 * @brief Reconstruct structure from torsion angles using constant bond info
 *
 * @param original_atoms std::vector<AtomCoordinate>
 * @param torsion_angles std::vector<float>
 * @return std::vector<AtomCoordinate>
 */
std::vector<AtomCoordinate> Nerf::reconstructWithConstants(
    std::vector<AtomCoordinate> original_atoms,
    std::vector<float> torsion_angles
) {
    std::vector<AtomCoordinate> reconstructed_atoms;
    // save three first atoms
    reconstructed_atoms.push_back(original_atoms[0]);
    reconstructed_atoms.push_back(original_atoms[1]);
    reconstructed_atoms.push_back(original_atoms[2]);

    int total = original_atoms.size();
    std::vector<AtomCoordinate> prev_atoms;
    std::vector< std::vector<float> > prev_coord;
    AtomCoordinate curr_atom;
    std::string curr_bond_name;
    std::vector<float> curr_coord;

    for (int i = 0; i < (total-3); i++) {
        curr_atom = original_atoms[i+3];
        prev_atoms = {reconstructed_atoms[i], reconstructed_atoms[i+1], reconstructed_atoms[i+2]};
        prev_coord = extractCoordinates(prev_atoms);
        curr_bond_name = original_atoms[i+2].atom + "_TO_" + curr_atom.atom;
        curr_coord = place_atom(
            prev_coord, this->bond_lengths[curr_bond_name],
            this->bond_angles[curr_bond_name], torsion_angles[i]
        );
        AtomCoordinate reconstructed_atom = AtomCoordinate(
            curr_atom.atom, curr_atom.residue, curr_atom.chain,
            curr_atom.atom_index, curr_atom.residue_index,
            curr_coord
        );
        reconstructed_atoms.push_back(reconstructed_atom);
    }

    return reconstructed_atoms;
}

/**
 * @brief Reconstruct structure from torsion angles using amino-acid-specific bond info
 *
 * @param original_atoms std::vector<AtomCoordinate>
 * @param torsion_angles std::vector<float>
 * @param aa_bond_lengths
 * @param aa_bond_angles
 * @return std::vector<AtomCoordinate>
 */
std::vector<AtomCoordinate> Nerf::reconstructWithAASpecificAngles(
    std::vector<AtomCoordinate> original_atoms,
    std::vector<float> torsion_angles,
    std::map<std::string, float> aa_bond_lengths,
    std::map< std::string, std::map<std::string, float> > aa_bond_angles
) {

    std::vector<AtomCoordinate> reconstructed_atoms;

    // save three first atoms
    reconstructed_atoms.push_back(original_atoms[0]);
    reconstructed_atoms.push_back(original_atoms[1]);
    reconstructed_atoms.push_back(original_atoms[2]);

    int total = original_atoms.size();
    std::vector<AtomCoordinate> prev_atoms;
    std::vector<std::vector<float> > prev_coord;
    AtomCoordinate curr_atom;
    std::string curr_bond_name;
    float curr_bond_length;
    float curr_bond_angle;
    std::vector<float> curr_coord;

    for (int i = 0; i < (total - 3); i++) {
        curr_atom = original_atoms[i + 3];
        prev_atoms = {reconstructed_atoms[i], reconstructed_atoms[i + 1],
                      reconstructed_atoms[i + 2]};
        prev_coord = extractCoordinates(prev_atoms);
        curr_bond_name = original_atoms[i + 2].atom + "_TO_" + curr_atom.atom;
        curr_bond_length = aa_bond_lengths[curr_bond_name];
        curr_bond_angle = aa_bond_angles[curr_atom.residue][curr_bond_name];
        curr_coord = place_atom(
            prev_coord, curr_bond_length, curr_bond_angle, torsion_angles[i]
        );
        AtomCoordinate reconstructed_atom = AtomCoordinate(
            curr_atom.atom, curr_atom.residue, curr_atom.chain,
            curr_atom.atom_index, curr_atom.residue_index, curr_coord);
        reconstructed_atoms.push_back(reconstructed_atom);
    }
    return reconstructed_atoms;
}

std::vector<AtomCoordinate> Nerf::reconstructAminoAcid(
    std::vector<AtomCoordinate> original_atoms,
    std::vector<float> torsion_angles,
    AminoAcid& aa
) {
    // Declare output
    std::vector<AtomCoordinate> reconstructed_atoms;

    // save three first atoms
    reconstructed_atoms.push_back(original_atoms[0]);
    reconstructed_atoms.push_back(original_atoms[1]);
    reconstructed_atoms.push_back(original_atoms[2]);

    std::map<std::string, float> aa_bond_lengths = aa.bondLengths;
    std::map<std::string, float> aa_bond_angles = aa.bondAngles;
    std::map<std::string, std::vector<std::string> > aa_side_chain = aa.sideChain;

    int total = aa.atoms.size();
    std::vector<AtomCoordinate> prev_atoms(3);
    std::vector<std::vector<float> > prev_coord;
    AtomCoordinate curr_atom = AtomCoordinate(
        "", original_atoms[0].residue, original_atoms[0].chain,
        original_atoms[2].atom_index + 1, original_atoms[0].residue_index,
        std::vector<float>(3, 0.0)
    );
    std::string curr_bond_name, curr_angle_name;
    float curr_bond_length, curr_bond_angle;

    std::vector<float> curr_coord;

    for (int i = 0; i < (total - 3); i++) {
        // Get current atom's info
        curr_atom.atom_index = reconstructed_atoms[i + 2].atom_index + 1;
        curr_atom.atom = aa.atoms[i + 3];
        // Fill prev_atoms
        for (int j = 0; j < 3; j++) {
            prev_atoms[j] = findFirstAtom(reconstructed_atoms, aa_side_chain[curr_atom.atom][j]);
        }
        prev_coord = extractCoordinates(prev_atoms);
        // Bond name & Angle name
        curr_bond_name = prev_atoms[2].atom + "_" + curr_atom.atom;
        curr_angle_name = prev_atoms[1].atom + "_" + prev_atoms[2].atom + "_" + curr_atom.atom;

        curr_bond_length = aa.bondLengths[curr_bond_name];
        curr_bond_angle = aa.bondAngles[curr_angle_name];

        curr_coord = place_atom(
            prev_coord, curr_bond_length, curr_bond_angle, torsion_angles[i]
        );

        AtomCoordinate reconstructed_atom = AtomCoordinate(
            curr_atom.atom, curr_atom.residue, curr_atom.chain,
            curr_atom.atom_index, curr_atom.residue_index, curr_coord
        );
        reconstructed_atoms.push_back(reconstructed_atom);
    }
    return reconstructed_atoms;
}

std::vector<AtomCoordinate> Nerf::reconstructWithAAMaps(
    std::vector<AtomCoordinate> original_atoms,
    std::map< std::string, std::vector<std::string> > prev_atom_map,
    std::map<std::string, float> aa_torsion_angles,
    std::map<std::string, float> aa_bond_lengths,
    std::map<std::string, float> aa_bond_angles
) {
    std::map< std::string, AtomCoordinate > original_atom_map;
    for (auto atm : original_atoms) {
        original_atom_map[atm.atom] = atm;
    }

    std::vector<AtomCoordinate> reconstructed_atoms;

    // save three first atoms
    reconstructed_atoms.push_back(original_atoms[0]);
    reconstructed_atoms.push_back(original_atoms[1]);
    reconstructed_atoms.push_back(original_atoms[2]);

    int total = original_atoms.size();
    std::vector<std::string> prev_atom_names;


    std::vector<std::vector<float> > prev_coord;
    AtomCoordinate curr_atom;
    std::string curr_bond_length_key;
    std::string curr_bond_angle_key;
    std::string curr_torsion_angle_key;

    float curr_bond_length;
    float curr_bond_angle;
    float curr_torsion_angle;

    std::vector<float> curr_coord;

    for (int i = 0; i < (total - 3); i++) {
        curr_atom = original_atoms[i + 3];
        prev_atom_names = prev_atom_map[curr_atom.atom];
        std::vector<AtomCoordinate> prev_atoms;
        for (auto atm_name : prev_atom_names) {
            prev_atoms.push_back(original_atom_map[atm_name]);
        }
        prev_coord = extractCoordinates(prev_atoms);

        curr_bond_length_key = prev_atoms[2].atom + "_" + curr_atom.atom;
        curr_bond_angle_key = prev_atoms[1].atom + "_" + prev_atoms[2].atom + "_" + curr_atom.atom;
        curr_torsion_angle_key = prev_atoms[0].atom + "_" + prev_atoms[1].atom +
            "_" + prev_atoms[2].atom + "_" + curr_atom.atom;

        curr_bond_length = aa_bond_lengths[curr_bond_length_key];
        curr_bond_angle = aa_bond_angles[curr_bond_angle_key];
        // Found error here! 2021-10-15 01:03:03
        curr_torsion_angle = aa_torsion_angles[curr_torsion_angle_key];
        curr_coord = place_atom(
            prev_coord, curr_bond_length, curr_bond_angle, curr_torsion_angle
        );
        AtomCoordinate reconstructed_atom = AtomCoordinate(
            curr_atom.atom, curr_atom.residue, curr_atom.chain,
            curr_atom.atom_index, curr_atom.residue_index, curr_coord);
        reconstructed_atoms.push_back(reconstructed_atom);
    }
    return reconstructed_atoms;
}




/**
 * @brief Reconstruct structure from torsion angles using atom-specific bond info
 *
 * @param original_atoms std::vector<AtomCoordinate>
 * @param torsion_angles
 * @param atom_bond_angles
 * @return std::vector<AtomCoordinate>
 */
std::vector<AtomCoordinate> Nerf::reconstructWithDynamicAngles(
    std::vector<AtomCoordinate> original_atoms,
    std::vector<float> torsion_angles,
    std::vector<float> atom_bond_lengths,
    std::vector<float> atom_bond_angles
) {
    std::vector<AtomCoordinate> reconstructed_atoms;

    // save three first atoms
    reconstructed_atoms.push_back(original_atoms[0]);
    reconstructed_atoms.push_back(original_atoms[1]);
    reconstructed_atoms.push_back(original_atoms[2]);

    int total = original_atoms.size();
    std::vector<AtomCoordinate> prev_atoms;
    std::vector<std::vector<float> > prev_coord;
    AtomCoordinate curr_atom;
    std::string curr_bond_name;
    float curr_bond_length;
    float curr_bond_angle;
    std::vector<float> curr_coord;

    for (int i = 0; i < (total - 3); i++) {
      curr_atom = original_atoms[i + 3];
      prev_atoms = {reconstructed_atoms[i], reconstructed_atoms[i + 1],
                    reconstructed_atoms[i + 2]};
      prev_coord = extractCoordinates(prev_atoms);
      curr_bond_name = original_atoms[i + 2].atom + "_TO_" + curr_atom.atom;
      curr_bond_length = atom_bond_lengths[i + 2];
      curr_bond_angle = atom_bond_angles[i + 1];
      curr_coord = place_atom(prev_coord, curr_bond_length, curr_bond_angle,
                              torsion_angles[i]);
      AtomCoordinate reconstructed_atom = AtomCoordinate(
          curr_atom.atom, curr_atom.residue, curr_atom.chain,
          curr_atom.atom_index, curr_atom.residue_index, curr_coord);
      reconstructed_atoms.push_back(reconstructed_atom);
    }
    return reconstructed_atoms;

}

std::vector<AtomCoordinate> Nerf::reconstructWithDynamicAngles(
    std::vector<AtomCoordinate> original_atoms,
    std::vector<float> torsion_angles,
    std::vector<float> atom_bond_angles
) {
    std::vector<AtomCoordinate> reconstructed_atoms;

    // save three first atoms
    reconstructed_atoms.push_back(original_atoms[0]);
    reconstructed_atoms.push_back(original_atoms[1]);
    reconstructed_atoms.push_back(original_atoms[2]);

    int total = original_atoms.size();
    std::vector<AtomCoordinate> prev_atoms;
    std::vector<std::vector<float> > prev_coord;
    AtomCoordinate curr_atom;
    std::string curr_bond_name;
    float curr_bond_length;
    float curr_bond_angle;
    std::vector<float> curr_coord;

    for (int i = 0; i < (total - 3); i++) {
      curr_atom = original_atoms[i + 3];
      prev_atoms = {reconstructed_atoms[i], reconstructed_atoms[i + 1],
                    reconstructed_atoms[i + 2]};
      prev_coord = extractCoordinates(prev_atoms);
      curr_bond_name = original_atoms[i + 2].atom + "_TO_" + curr_atom.atom;
      curr_bond_length = this->bond_lengths[curr_bond_name];
      curr_bond_angle = atom_bond_angles[i + 1];
      curr_coord = place_atom(prev_coord, curr_bond_length, curr_bond_angle,
                              torsion_angles[i]);
      AtomCoordinate reconstructed_atom = AtomCoordinate(
          curr_atom.atom, curr_atom.residue, curr_atom.chain,
          curr_atom.atom_index, curr_atom.residue_index, curr_coord);
      reconstructed_atoms.push_back(reconstructed_atom);
    }
    return reconstructed_atoms;
}


std::vector<AtomCoordinate> Nerf::reconstructWithBreaks(
    std::vector<AtomCoordinate> original_atoms,
    std::vector<float> torsion_angles,
    std::vector<float> atom_bond_angles,
    std::vector<int> break_indices
) {
    int total = original_atoms.size();
    std::vector<AtomCoordinate> reconstructed_atoms;
    int breakpoint;
    int next_breakpoint;
    std::vector<AtomCoordinate> prev_atoms;
    std::vector<std::vector<float> > prev_coord;
    AtomCoordinate curr_atom;
    std::string curr_bond_name;
    float curr_bond_length;
    float curr_bond_angle;
    std::vector<float> curr_coord;

    // save three first atoms
    for (size_t k = 0; k < break_indices.size(); k++) {
        breakpoint = break_indices[k];
        if (k != (break_indices.size() - 1)) {
            next_breakpoint = break_indices[k + 1];
        } else {
            next_breakpoint = total;
        }
        reconstructed_atoms.push_back(original_atoms[breakpoint]);
        reconstructed_atoms.push_back(original_atoms[breakpoint + 1]);
        reconstructed_atoms.push_back(original_atoms[breakpoint + 2]);

        for (int i = breakpoint; i < (next_breakpoint - 3); i++) {
            curr_atom = original_atoms[i + 3];
            prev_atoms = { reconstructed_atoms[i], reconstructed_atoms[i + 1],
                          reconstructed_atoms[i + 2] };
            prev_coord = extractCoordinates(prev_atoms);
            curr_bond_name = original_atoms[i + 2].atom + "_TO_" + curr_atom.atom;
            curr_bond_length = this->bond_lengths[curr_bond_name];
            curr_bond_angle = atom_bond_angles[i + 1];
            curr_coord = place_atom(prev_coord, curr_bond_length, curr_bond_angle,
                torsion_angles[i]);
            AtomCoordinate reconstructed_atom = AtomCoordinate(
                curr_atom.atom, curr_atom.residue, curr_atom.chain,
                curr_atom.atom_index, curr_atom.residue_index, curr_coord);
            reconstructed_atoms.push_back(reconstructed_atom);
        }
    }
    return reconstructed_atoms;
}

/**
 * @brief Reconstruct with relative positions
 *
 * @param original_atoms
 * @param relative_position
 * @return std::vector<AtomCoordinate>
 */
std::vector<AtomCoordinate> reconstrutWithRelativePositions(
    std::vector<AtomCoordinate> /* original_atoms */,
    std::vector<AtomCoordinate> /* relative_position */
){
    std::vector<AtomCoordinate> reconstructed_atoms;
    // Reconstruct with relative postion
    // Just adding them?

    return reconstructed_atoms;
}

// std::vector< std::vector<float> > Nerf::calculateCoordinates(
//     std::vector<AtomCoordinate> first_atoms,
//     std::vector<AtomCoordinate> last_atoms,
//     std::vector<float> torsion_angles,
//     std::vector<float> atom_bond_angles,
//     std::vector<float> atom_bond_lengths
// ) {
//     std::vector< std::vector<float> > output;

//     std::vector<std::vector<float> > prev_coord = extractCoordinates(first_atoms);
//     AtomCoordinate curr_atom;
//     std::string curr_bond_name;
//     float curr_bond_length;
//     float curr_bond_angle;
//     std::vector<float> curr_coord;

//     // save three first atoms
//     for (int i = 0; i < 3; i++) {
//         output.push_back(first_atoms[i]);
//     }

//     for (int i = 0; i < (total - 3); i++) {
//       curr_atom = atoms[i + 3];
//       prev_atoms = {reconstructed_atoms[i], reconstructed_atoms[i + 1],
//                     reconstructed_atoms[i + 2]};
//       prev_coord = extractCoordinates(prev_atoms);
//       curr_bond_name = atoms[i + 2].atom + "_TO_" + curr_atom.atom;
//       curr_bond_length = this->bond_lengths[curr_bond_name];
//       curr_bond_angle = atom_bond_angles[i + 1];
//       curr_coord = place_atom(prev_coord, curr_bond_length, curr_bond_angle,
//                               torsion_angles[i]);
//       AtomCoordinate reconstructed_atom = AtomCoordinate(
//           curr_atom.atom, curr_atom.residue, curr_atom.chain,
//           curr_atom.atom_index, curr_atom.residue_index, curr_coord);
//       reconstructed_atoms.push_back(reconstructed_atom);
//     }
//     return reconstructed_atoms;
// }
// 2022-04-14 21:33:03
// TODO: WRITE A METHOD THAT CAN RECONSTRUCT ONLY WITH 3 PREV ATOMS AND ANGLES


std::vector<AtomCoordinate> Nerf::reconstructWithReversed(
    std::vector<AtomCoordinate> original_atoms,
    std::vector<float> torsion_angles,
    std::vector<float> atom_bond_angles
) {
    std::vector<AtomCoordinate> reconstructed_atoms;

    std::reverse(original_atoms.begin(), original_atoms.end());
    std::reverse(torsion_angles.begin(), torsion_angles.end());
    std::reverse(atom_bond_angles.begin(), atom_bond_angles.end());

    // save three first atoms
    reconstructed_atoms.push_back(original_atoms[0]);
    reconstructed_atoms.push_back(original_atoms[1]);
    reconstructed_atoms.push_back(original_atoms[2]);

    int total = original_atoms.size();
    std::vector<AtomCoordinate> prev_atoms;
    std::vector<std::vector<float> > prev_coord;
    AtomCoordinate curr_atom;
    std::string curr_bond_name;
    float curr_bond_length;
    float curr_bond_angle;
    std::vector<float> curr_coord;

    for (int i = 0; i < (total - 3); i++) {
        curr_atom = original_atoms[i + 3];
        prev_atoms = { reconstructed_atoms[i], reconstructed_atoms[i + 1],
                      reconstructed_atoms[i + 2] };
        prev_coord = extractCoordinates(prev_atoms);
        curr_bond_name = curr_atom.atom  + "_TO_" +original_atoms[i + 2].atom;
        curr_bond_length = this->bond_lengths[curr_bond_name];
        curr_bond_angle = atom_bond_angles[i + 1];
        curr_coord = place_atom(prev_coord, curr_bond_length, curr_bond_angle,
            torsion_angles[i]);
        AtomCoordinate reconstructed_atom = AtomCoordinate(
            curr_atom.atom, curr_atom.residue, curr_atom.chain,
            curr_atom.atom_index, curr_atom.residue_index, curr_coord);
        reconstructed_atoms.push_back(reconstructed_atom);
    }

    // reverse the output
    std::reverse(reconstructed_atoms.begin(), reconstructed_atoms.end());

    return reconstructed_atoms;
}


/**
 * @brief Write coordinate and bond info into csv file
 *
 * @param coord_list std::vector<AtomCoordinate>
 * @param output_path std::string, path for output
 * @param sep std::string, separation character (default: comma)
 */
void Nerf::writeInfoForChecking(
    std::vector<AtomCoordinate> coord_list,
    std::string output_path, std::string sep
) {
    // setup file
    std::ofstream outfile;
    const int total = coord_list.size();
    // Things to write
    // CSV HEADER:
    //  AtomIndex,Atom,ResidueIndex,Residue,X,Y,Z,BondName,BondAngle,BondLength

    AtomCoordinate curr_atm;
    AtomCoordinate prev_atm, next_atm;
    std::string bond_name;
    std::string NA = "NA";
    float bond_angle, bond_length;
    std::string header = "AtomIndex,Atom,ResidueIndex,Residue,Chain,X,Y,Z,BondName,BondAngle,BondLength";

    outfile.open(output_path, std::ios::out);
    // Insert header
    outfile << header << "\n";

    for (int i = 0; i < total; i++) {
        //
        curr_atm = coord_list[i];

        if (i == 0) {
            outfile << curr_atm.atom_index << sep << curr_atm.atom << sep;
            outfile << curr_atm.residue_index << sep << curr_atm.residue << sep;
            outfile << curr_atm.chain << sep;
            outfile << curr_atm.coordinate[0] << sep << curr_atm.coordinate[1] << sep;
            outfile << curr_atm.coordinate[2] << sep;
            outfile << NA << sep << NA << sep << NA << "\n";
        } else {
            // WARNING:
            // FIXME: FIX THIS --> curr_atm to next atm
            prev_atm = coord_list[i - 1];
            bond_name = prev_atm.atom + "_TO_" + curr_atm.atom;
            bond_length = distance(prev_atm.coordinate, curr_atm.coordinate);
            // check error conditions
            // 01. same atom name with previous atom (ex: CA-CA bond)
            // --> SKIP PRINT
            if (prev_atm.atom == curr_atm.atom) {
                continue;
            }
            // print coordinate
            outfile << curr_atm.atom_index << sep << curr_atm.atom << sep;
            outfile << curr_atm.residue_index << sep << curr_atm.residue << sep;
            outfile << curr_atm.chain << sep;
            outfile << curr_atm.coordinate[0] << sep << curr_atm.coordinate[1] << sep;
            outfile << curr_atm.coordinate[2] << sep;
            // 02. different chain
            if (prev_atm.chain != curr_atm.chain) {
                outfile << NA << sep << NA << NA << "\n";
                continue;
            }

            if (i == (total-1)) {
                outfile << bond_name << sep << NA << sep << bond_length << "\n";
            } else {
                next_atm = coord_list[i + 1];
                // handle error condition 01.
                if (next_atm.atom == curr_atm.atom) {
                    next_atm = coord_list[i + 2];
                }
                bond_angle = angle(
                    prev_atm.coordinate, curr_atm.coordinate,
                    next_atm.coordinate
                );
                outfile << bond_name << sep << bond_angle << sep;
                outfile << bond_length << "\n";
            }
        }
    }
    outfile.close();
}


void Nerf::writeCoordinatesBinary(
    std::vector<AtomCoordinate> coord_list, std::string output_path
) {
    // setup file
    std::ofstream outfile;
    const int total = coord_list.size();
    // Things to write
    // CSV HEADER:
    //  AtomIndex,Atom,ResidueIndex,Residue,X,Y,Z,BondName,BondAngle,BondLength

    AtomCoordinate curr_atm;
    // std::string header = "AtomIndex,Atom,ResidueIndex,Residue,Chain,X,Y,Z";

    outfile.open(output_path, std::ios::out | std::ios::binary);

    for (int i = 0; i < total; i++) {
        curr_atm = coord_list[i];
        outfile.write((char*)&curr_atm.atom, 1);
        outfile.write((char*)&curr_atm.coordinate[0], 4);
        outfile.write((char*)&curr_atm.coordinate[1], 4);
        outfile.write((char*)&curr_atm.coordinate[2], 4);
    }
    outfile.close();
}



std::vector<float> Nerf::getBondAngles(std::vector<AtomCoordinate> original_atoms) {
    // define variables
    std::vector<float> output;
    AtomCoordinate curr_atm, prev_atm, next_atm;
    float bond_angle;
    const int total = original_atoms.size();
    for (int i = 1; i < (total - 1); i++) {
      //
        curr_atm = original_atoms[i];
        prev_atm = original_atoms[i - 1];
        next_atm = original_atoms[i + 1];
        bond_angle = angle(prev_atm.coordinate, curr_atm.coordinate,
                           next_atm.coordinate);
        output.push_back(bond_angle);
    }
    // ADD LAST ANGLE
    return output;
}

std::vector<float> Nerf::getBondLengths(std::vector<AtomCoordinate> original_atoms) {
    // define variables
    std::vector<float> output;
    AtomCoordinate curr_atm, next_atm;
    float bond_length;
    const int total = original_atoms.size();
    for (int i = 0; i < (total - 1); i++) {
      //
      curr_atm = original_atoms[i];
      next_atm = original_atoms[i + 1];
      bond_length = distance(curr_atm.coordinate, next_atm.coordinate);
      output.push_back(bond_length);
    }
    return output;
}

/**
 * @brief Identify breaks in chain and return the indices of break point
 *
 * @param original_atoms A vector of AtomCoordinate
 * @param cutoff Standard to define a break. Default value = 3
 * @return std::vector<int>
 */
std::vector<int> Nerf::identifyBreaks(std::vector<AtomCoordinate> original_atoms, float cutoff) {
    // define variables
    std::vector<int> output = {0};
    std::vector<float> bondLengths = this->getBondLengths(original_atoms);

    // No need to start with 0
    for (size_t i = 1; i < bondLengths.size(); i++) {
        // If length is bigger than the cutoff, add breakpoint
        if (bondLengths[i] > cutoff) {
            // bondLenght[i] = dist between original_atoms[i], [i+1]
            // append i+1 to the output
            if (i <= (original_atoms.size() - 3)) {
                output.push_back(i + 1);
                // std::cout << i + 1 << std::endl;
            }
        }
    }
    return output;
}

// std::vector<float> compute_positions(
//     std::vector<AtomCoordinate> coord_list
// ) {

//     std::vector<float> output {0.0, 0.0};
//     return output;
// }
