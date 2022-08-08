/**
 * File: pdb_reader.cpp
 * Project: src
 * Created: 2021-01-04 17:31:03
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     Functions and a class for reading PDB files
 *     Deprecated. Use StructureReader in "input.h" instead.
 * ---
 * Last Modified: 2022-08-04 16:56:54
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */

#include "pdb_reader.h"


/**
 * Read title
 * @brief Read title of pdb file
 * @return std::string
 */
std::string PDBReader::readTitle() {
    std::string line;
    std::string rec_type;
    std::string empty_title = "";
    bool title_found = false;
    this->infile.clear();
    this->infile.seekg(0, std::ios::beg);
    while (std::getline(this->infile, line)) {
        rec_type = line.substr(0, 6);
        if (rec_type == "TITLE ") {
            if (!title_found) {
                title_found = true;
                this->title += line.substr(10);
            } else {
                this->title += line.substr(8);
            }
        }
    }
    if (this->title.empty()) {
        return empty_title;
    } else {
        return this->title;
    }
}


std::string remove_whitespaces(std::string input) {
    input.erase(
        std::remove_if(input.begin(), input.end(), isspace), input.end()
    );
    return input;
}

//
float* PDBReader::getOneAtomCoordinate(int serial_number) {
    static float output[3];
    std::string line;
    // variables for parsing line
    std::string rec_type;
    std::string str_serial_num;
    std::string str_x, str_y, str_z;
    int curr_serial_num;
    // This works fine 2021-01-06 12:27:59
    // read line by line
    while (std::getline(this->infile, line)) {
        rec_type = line.substr(0, 6);
        if (rec_type == "ATOM  ") {
            // NOTE: index for each item
            // serial number: 6-11 --> length: 5
            // coordinate: x(30-38), y(38-46), z(46-54) --> length: 8
            // str.substr(pos, len)
            str_serial_num = line.substr(6, 5);
            str_serial_num = remove_whitespaces(str_serial_num);
            curr_serial_num = std::stoi(str_serial_num);
            if (curr_serial_num == serial_number) {
                str_x = line.substr(30, 8);
                str_y = line.substr(38, 8);
                str_z = line.substr(46, 8);
                // delete whitespaces
                str_x = remove_whitespaces(str_x);
                str_y = remove_whitespaces(str_y);
                str_z = remove_whitespaces(str_z);
                // parse into float
                output[0] = std::stof(str_x);
                output[1] = std::stof(str_y);
                output[2] = std::stof(str_z);
                return output;
            }
        } else {
            continue;
        }
    }
    return output;
}


std::vector< std::vector<float> > PDBReader::getAllAtomCoordinates(
    std::string atom_name
) {
    std::vector< std::vector<float> > output;
    // variables for parsing line
    std::string line;
    std::string rec_type;
    std::string curr_atom_name;
    std::string str_x, str_y, str_z;
    float x, y, z;

    // FIXME: Code duplicate. Need to be fixed!
    while (std::getline(this->infile, line)) {
        rec_type = line.substr(0, 6);
        if (rec_type == "ATOM  ") {
            curr_atom_name = line.substr(12, 4);
            curr_atom_name = remove_whitespaces(curr_atom_name);
            if (curr_atom_name == atom_name) {
                str_x = line.substr(30, 8);
                str_y = line.substr(38, 8);
                str_z = line.substr(46, 8);
                // delete whitespaces
                str_x = remove_whitespaces(str_x);
                str_y = remove_whitespaces(str_y);
                str_z = remove_whitespaces(str_z);
                // parse into float
                x = std::stof(str_x);
                y = std::stof(str_y);
                z = std::stof(str_z);
                std::vector<float> curr_coord { x, y, z };
                output.push_back(curr_coord);
            }
        } else {
            continue;
        }
    }
    return output;
}

std::vector< std::vector<float> > PDBReader::getAllAtomCoordinates(
    std::vector<std::string> atom_list
) {
    std::vector< std::vector<float> > output;
    // variables for parsing line
    std::string line;
    std::string rec_type;
    std::string curr_atom_name;
    std::string str_x, str_y, str_z;
    float x, y, z;
    std::vector<std::string>::iterator str_it;

    // FIXME: Code duplicate. Need to be fixed!
    while (std::getline(this->infile, line)) {
        rec_type = line.substr(0, 6);
        if (rec_type == "ATOM  ") {
            curr_atom_name = line.substr(12, 4);
            curr_atom_name = remove_whitespaces(curr_atom_name);
            // check curr_atom_name in atom_list
            str_it = std::find(atom_list.begin(), atom_list.end(), curr_atom_name);
            if (str_it != atom_list.end()) {
                str_x = line.substr(30, 8);
                str_y = line.substr(38, 8);
                str_z = line.substr(46, 8);
                // delete whitespaces
                str_x = remove_whitespaces(str_x);
                str_y = remove_whitespaces(str_y);
                str_z = remove_whitespaces(str_z);
                // parse into float
                x = std::stof(str_x);
                y = std::stof(str_y);
                z = std::stof(str_z);
                std::vector<float> curr_coord{ x, y, z };
                output.push_back(curr_coord);
            }
        }
        else {
            continue;
        }
    }
    return output;
}

void PDBReader::writeCoordinates(std::string atom_name, std::string file_path) {
    // new ofstream
    std::ofstream output(
        file_path,
        std::ios_base::binary|std::ios_base::out
    );
    std::vector<std::vector<float>> coordinates;
    // To write in binary file
    const char* cstr_atom = atom_name.c_str();

    coordinates = this->getAllAtomCoordinates(atom_name);
    for (std::vector<float> curr_coord : coordinates) {
        output.write(cstr_atom, 2); // 2 char: CA -> 2 bytes
        for (int i = 0; i < 3; i++) { // x, y, z
            if (output.good()) {
                output.write((char*) &curr_coord[i], sizeof(float)); // 4 bytes
            }
        }
    }
    output.close();
}

float _extractFloatFromPDBLine(std::string line, int start, int len) {
    std::string str_x = line.substr(start, len);
    str_x = remove_whitespaces(str_x);
    return std::stof(str_x);
}

int _extractIntFromPDBLine(std::string line, int start, int len) {
    std::string str_x = line.substr(start, len);
    str_x = remove_whitespaces(str_x);
    return std::stoi(str_x);
}

/**
 * @brief Load an atom coordinate with given line
 *
 * @param line A string of line in pdb file
 * @return AtomCoordinate
 */
AtomCoordinate PDBReader::loadOneAtomCoordinate(std::string line) {
    // variables for parsing line
    std::string curr_atom_name;
    std::string curr_residue;
    std::string curr_chain_id;
    int atom_index, residue_index;
    float x, y, z;
    float tf;

    curr_atom_name = line.substr(12, 4);
    curr_atom_name = remove_whitespaces(curr_atom_name);
    curr_residue = line.substr(17, 3);
    curr_residue = remove_whitespaces(curr_residue);
    curr_chain_id = line.substr(21, 1);

    // parse into integer
    atom_index = _extractIntFromPDBLine(line, 6, 5);
    residue_index = _extractIntFromPDBLine(line, 22, 4);
    // parse into float
    x = _extractFloatFromPDBLine(line, 30, 8);
    y = _extractFloatFromPDBLine(line, 38, 8);
    z = _extractFloatFromPDBLine(line, 46, 8);
    tf = _extractFloatFromPDBLine(line, 60, 6);

    AtomCoordinate output = AtomCoordinate(
        curr_atom_name, curr_residue, curr_chain_id,
        atom_index, residue_index, x, y, z
    );

    output.setTempFactor(tf);

    return output;

};

/**
 * @brief Load coordinates with atom list
 *
 * @param atom_list
 * @return std::vector<AtomCoordinate>
 */
std::vector<AtomCoordinate> PDBReader::loadAllAtomCoordinates(
    std::vector<std::string> atom_list
) {
    // declare variables for loading
    std::string line;
    std::string rec_type;
    std::string curr_atom_name;
    std::vector<std::string>::iterator str_it;
    AtomCoordinate temp_atom_coord;
    std::vector<AtomCoordinate> output;

    // Iterate line by line
    while (std::getline(this->infile, line)) {
        rec_type = line.substr(0, 6);
        if (rec_type == "ATOM  ") {
            curr_atom_name = line.substr(12, 4);
            curr_atom_name = remove_whitespaces(curr_atom_name);
            // check curr_atom_name in atom_list
            str_it = std::find(atom_list.begin(), atom_list.end(), curr_atom_name);
            if (str_it != atom_list.end()) {
                temp_atom_coord = loadOneAtomCoordinate(line);
                output.push_back(temp_atom_coord);
            }
        }
        else {
            continue;
        }
    }
    return output;
}

/**
 * @brief Load all atoms including the ones from side chain
 *
 * @return std::vector<AtomCoordinate>
 */
std::vector<AtomCoordinate> PDBReader::loadAllAtomCoordinatesWithSideChain() {
    // declare variables for loading
    std::string line;
    std::string rec_type;
    AtomCoordinate temp_atom_coord;
    std::vector<AtomCoordinate> output;

    // Iterate line by line
    while (std::getline(this->infile, line)) {
        rec_type = line.substr(0, 6);
        if (rec_type == "ATOM  ") {
            temp_atom_coord = loadOneAtomCoordinate(line);
            output.push_back(temp_atom_coord);
        }
        else {
            continue;
        }
    }
    return output;
}
