/**
 * File: input.h
 * Project: src
 * Created: 2022-08-04 15:02:17
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     This code is written as part of project "src".
 * References:
 *     GemmiWrapper.cpp in https://github.com/steineggerlab/foldseek
 * ---
 * Last Modified: 2022-08-04 16:42:32
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2022 Hyunbin Kim, All rights reserved
 */

#pragma once
#include "atom_coordinate.h"
#include <stdexcept>
// Gemmi
#include "gemmi/mmread.hpp"
#include "gemmi/input.hpp"
#include "gemmi/gz.hpp"

class StructureReader {
private:
    void updateStructure(void* void_st, std::string& filename);
public:
    StructureReader(/* args */){};
    ~StructureReader(){};
    std::string filepath;
    std::string title;
    std::vector<AtomCoordinate> atoms;
    bool loadFromBuffer(const char* buffer, size_t bufferSize, std::string& name);
    bool load(std::string& filename);
    bool readBackboneAtoms(std::vector<AtomCoordinate>& backboneAtoms);
    bool readAllAtoms(std::vector<AtomCoordinate>& allAtoms);
};
