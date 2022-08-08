/**
 * File: atom_coordinate.cpp
 * Project: src
 * Created: 2021-01-18 12:53:34
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     The data type to handle atom coordinate comes here.
 * ---
 * Last Modified: 2022-08-08 21:12:28
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */

#include "atom_coordinate.h"

// constructor
AtomCoordinate::AtomCoordinate(){}

/**
 * @brief Construct a new Atom Coordinate:: Atom Coordinate object
 *
 * @param a A string for atom name
 * @param r A string for residue name
 * @param ai An integer for atom index
 * @param ri An integer for residue index
 * @param x A float for x coordinate
 * @param y A float for y coordinate
 * @param z A float for z coordinate
 */
AtomCoordinate::AtomCoordinate(
    std::string a, std::string r, std::string c,
    int ai, int ri, float x, float y, float z
): atom(a), residue(r), chain(c), atom_index(ai), residue_index(ri) {
    this->coordinate = {x, y, z};
    this->check3dCoordinate();
}

/**
 * @brief Construct a new Atom Coordinate:: Atom Coordinate object
 *
 * @param a A string for atom name
 * @param r A string for residue name
 * @param ai An integer for atom index
 * @param ri An integer for residue index
 * @param coord A float vector for x,y,z coordinates.
 */
AtomCoordinate::AtomCoordinate(
    std::string a, std::string r, std::string c,
    int ai, int ri, std::vector<float> coord
): atom(a), residue(r), chain(c), atom_index(ai), residue_index(ri), coordinate(coord) {
    this->check3dCoordinate();
};


AtomCoordinate::~AtomCoordinate(){}

/**
 * @brief Check if 3 coordinates given.
 *
 */
void AtomCoordinate::check3dCoordinate() {
    if (this->coordinate.size() != 3) {
        throw(std::out_of_range(
            "Unexpected length of atom coordinates: 3 floats should be given"
        ));
    }
}

bool AtomCoordinate::operator==(const AtomCoordinate& other) const {
    return (
        (this->atom == other.atom) &&
        (this->atom_index == other.atom_index) &&
        (this->residue == other.residue) &&
        (this->residue_index == other.residue_index) &&
        (this->chain == other.chain) &&
        (this->coordinate == other.coordinate)
    );
}
bool AtomCoordinate::operator!=(const AtomCoordinate& other) const {
    return !(*this == other);
}

bool AtomCoordinate::isBackbone() {
    return ((this->atom == "N") ||(this->atom == "CA") ||(this->atom == "C"));
}

void AtomCoordinate::print(int option) {
    std::cout << "Atom: " << this->atom << std::endl;
    if (option != 0) {
        std::cout << "Residue: " << this->residue << std::endl;
        std::cout << "Chain: " << this->chain << std::endl;
        std::cout << "Atom Index: " << this->atom_index << std::endl;
        std::cout << "Residue Index: " << this->residue_index << std::endl;
        if (option == 2) {
            std::cout << "Coordinate: ";
            for (int i = 0; i < this->coordinate.size(); i++) {
                std::cout << this->coordinate[i] << " ";
            }
            std::cout << std::endl;
        }
    }
}

/**
 * @brief Extracts coordinates from AtomCoordinate vector
 *
 * @param atoms A vector of AtomCoordinate
 * @return std::vector< std::vector<float> >
 */
std::vector< std::vector<float> > extractCoordinates(
    std::vector<AtomCoordinate>& atoms
) {
    std::vector<std::vector<float>> output;
    for (AtomCoordinate curr_atm : atoms) {
        output.push_back(curr_atm.coordinate);
    }
    return output;
}

std::vector<AtomCoordinate> extractChain(
    std::vector<AtomCoordinate>& atoms, std::string chain
) {
    std::vector<AtomCoordinate> output;
    AtomCoordinate next_atm;
    AtomCoordinate curr_atm;
    int total = atoms.size();
    for (int i = 0; i < total; i++) {
        curr_atm = atoms[i];
        if (i < (total-1)) {
            next_atm = atoms[i + 1];
            if (next_atm.atom == curr_atm.atom) {
                continue;
            }
        }
        if (curr_atm.chain == chain) {
            output.push_back(curr_atm);
        }
    }
    return output;
}

void printAtomCoordinateVector(std::vector<AtomCoordinate>& atoms, int option) {
    for (AtomCoordinate curr_atm : atoms) {
        curr_atm.print(option);
    }
}

std::vector<AtomCoordinate> filterBackbone(std::vector<AtomCoordinate>& atoms) {
    std::vector<AtomCoordinate> output;
    for (AtomCoordinate curr_atm : atoms) {
        if (curr_atm.isBackbone()) {
            output.push_back(curr_atm);
        }
    }
    return output;
}

std::vector<AtomCoordinate> weightedAverage(
    std::vector<AtomCoordinate>& origAtoms, std::vector<AtomCoordinate>& revAtoms
) {
    std::vector<AtomCoordinate> output;
    AtomCoordinate curr_atm;
    AtomCoordinate rev_atm;
    int total = origAtoms.size();
    for (int i = 0; i < total; i++) {
        curr_atm = origAtoms[i];
        rev_atm = revAtoms[i];
        std::vector<float> coord = {
            ((curr_atm.coordinate[0] * (float)(total - i)) + (rev_atm.coordinate[0] * (float)i)) / (float)total,
            ((curr_atm.coordinate[1] * (float)(total - i)) + (rev_atm.coordinate[1] * (float)i)) / (float)total,
            ((curr_atm.coordinate[2] * (float)(total - i)) + (rev_atm.coordinate[2] * (float)i)) / (float)total
        };
        output.push_back(AtomCoordinate(
            curr_atm.atom, curr_atm.residue, curr_atm.chain,
            curr_atm.atom_index, curr_atm.residue_index, coord
        ));
    }
    return output;
}

// 2021-01-18 13:53:00 TESTED.

int writeAtomCoordinatesToPDB(
    std::vector<AtomCoordinate>& atoms, std::string title, std::string pdb_path
) {
    std::ofstream pdb_file;
    pdb_file.open(pdb_path);
    if (!pdb_file.is_open()) {
        std::cout << "Error: Cannot open file: " << pdb_path << std::endl;
        return 1;
    }
    // Write title
    // 2022-07-20 05:15:06 Currently NOT WORKING.
    if (title != "") {
        int title_line_num = (int)(ceil((title.length() - 70) / 72.0) + 1);
        // Split title into lines of 70 characters.
        std::vector<std::string> title_per_line;
        for (int i = 0; i < title_line_num; i++) {
            if (i == 0) {
                title_per_line.push_back(title.substr(0, 70));
            } else {
                title_per_line.push_back(title.substr(i * 72 - 2, 72));
            }
        }
        for (int i = 0; i < title_line_num; i++) {
            if (i == 0) {
                pdb_file << "TITLE     " << title_per_line[i] << "\n";
            } else {
                pdb_file << "TITLE   " << title_per_line[i] << std::endl;
            }
        }
    }

    int total = atoms.size();
    std::string residue;
    int residue_index;
    for (int i = 0; i < total; i++) {
        pdb_file << "ATOM  "; // 1-4 ATOM
        pdb_file << std::setw(5) << i + 1; // 7-11
        pdb_file << " "; // 12
        if (atoms[i].atom.size() == 4) {
            pdb_file << std::setw(4) << std::left << atoms[i].atom; // 13-16
        } else {
            pdb_file << " ";
            pdb_file << std::setw(3) << std::left << atoms[i].atom; // 13-16
        }
        pdb_file << " "; // 17
        pdb_file << std::setw(3) << std::right << atoms[i].residue; // 18-20
        pdb_file << " "; // 21
        pdb_file << atoms[i].chain; // 22
        pdb_file << std::setw(4) << atoms[i].residue_index; // 23-26
        pdb_file << "    "; // 27-30
        pdb_file << std::setw(8) << std::setprecision(3) << std::fixed << atoms[i].coordinate[0]; // 31-38
        pdb_file << std::setw(8) << std::setprecision(3) << std::fixed << atoms[i].coordinate[1]; // 39-46
        pdb_file << std::setw(8) << std::setprecision(3) << std::fixed << atoms[i].coordinate[2]; // 47-54
        pdb_file << "  1.00"; // 55-60
        pdb_file << std::setw(6) << std::setprecision(2) << std::fixed << atoms[i].tempFactor; // 61-66
        pdb_file << "          "; // 67-76
        // First one character from atom
        pdb_file << std::setw(2) << atoms[i].atom[0]; // 77-78
        pdb_file << "  \n"; // 79-80
        if (i == (total-1)) {
            // TER
            // 1-6 Record name "TER   "
            // 7-11 Atom serial number.
            // 18-20 Residue name.
            // 22 Chain identifier.
            // 23-26 Residue sequence number.
            pdb_file << "TER   " << std::setw(5) << total + 1 << "      ";
            pdb_file << std::setw(3) << std::right << atoms[i].residue;
            pdb_file << " " << atoms[i].chain;
            pdb_file << std::setw(4) << atoms[i].residue_index << std::endl;
        }
    }
    pdb_file.close();
    return 0;
}

std::vector< std::vector<AtomCoordinate> > splitAtomByResidue(
    const std::vector<AtomCoordinate>& atomCoordinates
) {
    std::vector< std::vector<AtomCoordinate> > output;
    std::vector<AtomCoordinate> currentResidue;

    for (int i = 0; i < atomCoordinates.size(); i++) {
        if (i == 0) {
            currentResidue.push_back(atomCoordinates[i]);
        } else if (i != (atomCoordinates.size() - 1)) {
            if (atomCoordinates[i].residue_index == atomCoordinates[i-1].residue_index) {
                currentResidue.push_back(atomCoordinates[i]);
            } else {
                output.push_back(currentResidue);
                currentResidue.clear();
                currentResidue.push_back(atomCoordinates[i]);
            }
        } else {
            currentResidue.push_back(atomCoordinates[i]);
            output.push_back(currentResidue);
        }
    }

    return output;
}

std::vector<std::string> getResidueNameVector(
    const std::vector<AtomCoordinate>& atomCoordinates
) {
    std::vector<std::string> output;
    // Unique residue names
    for (int i = 0; i < atomCoordinates.size(); i++) {
        if (i == 0) {
            output.push_back(atomCoordinates[i].residue);
        } else {
            if (atomCoordinates[i].residue_index != atomCoordinates[i-1].residue_index) {
                output.push_back(atomCoordinates[i].residue);
            }
        }
    }
    return output;
}

AtomCoordinate findFirstAtom(std::vector<AtomCoordinate>& atoms, std::string atom_name) {
    for (AtomCoordinate curr_atm : atoms) {
        if (curr_atm.atom == atom_name) {
            return curr_atm;
        }
    }
    return AtomCoordinate();
}

void setAtomIndexSequentially(std::vector<AtomCoordinate>& atoms, int start) {
    for (int i = 0; i < atoms.size(); i++) {
        atoms[i].atom_index = start + i;
    }
}

void removeAlternativePosition(std::vector<AtomCoordinate>& atoms) {
    // If there is an alternative position, remove it
    for (int i = 1; i < atoms.size(); i++) {
        if (atoms[i].atom == atoms[i-1].atom) {
            atoms.erase(atoms.begin() + i);
            i--;
        }
    }
}

std::vector<AtomCoordinate> getAtomsWithResidueIndex(
    std::vector<AtomCoordinate>& atoms, int residue_index,
    std::vector<std::string> atomNames
) {
    std::vector<AtomCoordinate> output;
    for (AtomCoordinate curr_atm : atoms) {
        if (curr_atm.residue_index == residue_index) {
            for (std::string atom_name : atomNames) {
                if (curr_atm.atom == atom_name) {
                    output.push_back(curr_atm);
                }
            }
        }
    }
    return output;
}

std::vector< std::vector<AtomCoordinate> > getAtomsWithResidueIndex(
    std::vector<AtomCoordinate>& atoms, std::vector<int> residue_index,
    std::vector<std::string> atomNames
) {
    std::vector< std::vector<AtomCoordinate> > output;
    for (int curr_index : residue_index) {
        std::vector<AtomCoordinate> curr_atoms = getAtomsWithResidueIndex(atoms, curr_index, atomNames);
        output.push_back(curr_atoms);
    }
    return output;
}