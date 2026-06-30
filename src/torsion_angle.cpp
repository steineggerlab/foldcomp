/**
 * File: torsion_angle.cpp
 * Project: foldcomp
 * Created: 2021-01-13 10:34:34
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     This code calculates torsion angles from the atomic coordinates.
 * Reference:
 *     1) "torsion.xyz.R" in Bio3D R package
 *     https://rdrr.io/cran/bio3d/man/torsion.xyz.html
 *     2) "XYZ.h" in pdbtools
 *     https://github.com/realbigws/PDB_Tool
 * ---
 * Last Modified: 2022-09-13 15:15:46
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
#include "torsion_angle.h"

#include "atom_coordinate.h"

#include <cstddef>
#include <cstdio>
#include <iostream>
#include <string>

 // Temp function
void print3DFloatVec(std::string name, const std::vector<float>& input) {
    std::cout << name << ": " << input[0] << ", " << input[1] << ", " << input[2] << std::endl;
}

TorsionTrace computeTorsionTraceRuntime(
    const float3d& atm1,
    const float3d& atm2,
    const float3d& atm3,
    const float3d& atm4
) {
    return computeTorsionTraceFromCoords(atm1, atm2, atm3, atm4);
}

std::vector<float> getTorsionFromXYZ(
    const std::vector<AtomCoordinate>& coordinates, int atm_inc
) {
    std::vector<float3d> coord = extractCoordinates(coordinates);
    return getTorsionFromXYZ(coord, atm_inc);
}

/**
 * @brief Get the torsion from xyz object
 *
 * @param coordinates a vector of float vectors
 * @param atm_inc an integer
 * @return std::vector<std::vector<float>>
 */
std::vector<float> getTorsionFromXYZ(
    const std::vector<float3d>& coordinates, int atm_inc = 1
) {
    std::vector<float> torsion_vector;
    if (coordinates.size() < 4 || atm_inc <= 0) {
        return torsion_vector;
    }
    for (size_t i = 0; i < (coordinates.size() - 3); i += atm_inc) {
        torsion_vector.push_back(computeTorsionAngleFromCoords(
            coordinates[i + 0],
            coordinates[i + 1],
            coordinates[i + 2],
            coordinates[i + 3]
        ));
    }
    return torsion_vector;
}


/**
 * WARNING: TEMPORARY FUNCTION
 * @brief Fill an array of double with float vector
 *
 * @param fv
 * @param output
 */
void float3dVectorToDoubleArray(const std::vector<float>& fv, double output[3]) {
    for (int i = 0; i < 3; i++) {
        output[i] = static_cast<double>(fv[i]);
    }
}


/**
 * @brief Save torsion angles to the output
 *
 * @param output_path
 * @param torsion
 */
void writeTorsionAngles(const std::string& file_path, const std::vector<float>& torsion) {
    FILE* output = fopen(file_path.c_str(), "w");
    if (output == NULL) {
        return;
    }
    for (float angle : torsion) {
        fprintf(output, "%f\n", angle);
    }
    fclose(output);
}

std::vector<float> readTorsionAngles(const std::string& file_path) {
    std::vector<float> output;
    FILE* input = fopen(file_path.c_str(), "r");
    if (input == NULL) {
        return output;
    }
    char line[128];
    while (fgets(line, sizeof(line), input) != NULL) {
        output.push_back(strtof(line, NULL));
    }
    fclose(input);
    return output;
}

// Encode backbone by encoding dihedrals by 2B each into intervals of 2 pi / 65536
// Current, the degree is not in radian and 360 will be used instead of 2 pi.

/**
 * @brief Encode float-based torsion angles to short
 *
 * @param torsions A float vector of torsion angles
 * @param n_bits Number of bits for encoding
 * @return std::vector<short>
 */
std::vector<short> encodeTorsionAnglesToShort(
    const std::vector<float>& torsions, unsigned int n_bits
) {
    std::vector<short> output;
    short s_angle;
    for (float f_angle: torsions) {
        s_angle = (short)(f_angle * pow(2, n_bits) / 360);
        output.push_back(s_angle);
    }
    return output;
}

/**
 * @brief Decode short-encoded torsion angles to float
 *
 * @param encoded_torsions A short vector of torsion angles
 * @param n_bits Number of bits for encoding
 * @return std::vector<float>
 */
std::vector<float> decodeEncodedTorsionAngles(
    const std::vector<short>& encoded_torsions, unsigned int n_bits
) {
    std::vector<float> output;
    float f_angle;
    for (short s_angle: encoded_torsions) {
        f_angle = (float)(s_angle);
        f_angle = f_angle * 360.0 / pow(2, n_bits);
        output.push_back(f_angle);
    }
    return output;
}
