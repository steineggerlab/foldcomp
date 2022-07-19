/**
 * File: pdb_reader.h
 * Project: foldcomp
 * Created: 2021-01-04 16:32:18
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     Functions and a class for reading PDB files
 * ---
 * Last Modified: 2022-07-20 01:53:38
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
#pragma once
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include "atom_coordinate.h"

class PDBReader {
private:
    /* data */

    bool is_loaded;
    void check_loaded();

public:
    /* properties */

    int pdb_line_length;
    std::string file_path;
    std::string delimiter;
    std::ifstream infile;
    std::string title;

    /* constructor & destructor */

    PDBReader():is_loaded(0), pdb_line_length(80), delimiter(" "), file_path("") {};
    PDBReader(std::string file_path):is_loaded(0), pdb_line_length(80), delimiter(" "), file_path("") {
        this->file_path = file_path;
    };
    ~PDBReader(){};

    /* methods */

    void load() {
        if (this->file_path != "") {
            this->infile.open(this->file_path);
        }
    };
    void load(std::string file_path) {
        this->infile.open(file_path);
    };

    std::string readTitle();
    void getTitle();
    float* getOneAtomCoordinate(int serial_number);
    std::vector<std::vector <float> > getAllAtomCoordinates(std::string atom_name);
    std::vector< std::vector<float> > getAllAtomCoordinates(
        std::vector<std::string> atom_list = std::vector<std::string>{"N","CA","C"}
    );
    // Added 2021-01-18 16:47:33
    // Support for AtomCoordinate class
    AtomCoordinate loadOneAtomCoordinate(std::string line);
    std::vector<AtomCoordinate> loadAllAtomCoordinates(
        std::vector<std::string> atom_list = std::vector<std::string>{"N","CA","C"}
    );

    // 2021-08-18 13:27:57
    std::vector<AtomCoordinate> loadAllAtomCoordinatesWithSideChain();

    void writeCoordinates(std::string atom_name, std::string file_path);
    // TODO: NEED TO IMPLEMENT THIS WITH UPDATED FORMAT FOR COORDINATES
    // void writeCoordinates(std::vector<std::string> atom_list, std::string file_path);
};
