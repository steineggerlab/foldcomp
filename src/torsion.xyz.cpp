/**
 * File: torsion.xyz.cpp
 * Project: src
 * Created: 2021-01-13 10:34:34
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     This code is written as part of project "src".
 *     This code calculates torsion angles from the atomic coordinates.
 * Reference:
 *     1) "torsion.xyz.R" in Bio3D R package
 *     https://rdrr.io/cran/bio3d/man/torsion.xyz.html
 *     2) "XYZ.h" in pdbtools
 *     https://github.com/realbigws/PDB_Tool
 * ---
 * Last Modified: 2022-06-02 19:31:43
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */

#include "torsion.xyz.h"



 // Temp function
void print3DFloatVec(std::string name, std::vector<float> input) {
    std::cout << name << ": " << input[0] << ", " << input[1] << ", " << input[2] << std::endl;
}

/**
 * @brief Get the torsion from xyz object
 *
 * @param coordinates a vector of float vectors
 * @param atm_inc an integer
 * @return std::vector<std::vector<float>>
 */
std::vector<float> getTorsionFromXYZForDebugging(
    std::vector< std::vector<float> > coordinates, int atm_inc = 1
) {
    std::vector<float> torsion_vector;
    std::vector<int> atm_inds{ 0,1,2,3 };
    while (atm_inds[3] < coordinates.size()) {
        // 00. Get 4 atom coordinates
        std::vector<float> atm_1 = coordinates[atm_inds[0]];
        std::vector<float> atm_2 = coordinates[atm_inds[1]];
        std::vector<float> atm_3 = coordinates[atm_inds[2]];
        std::vector<float> atm_4 = coordinates[atm_inds[3]];
        // 01. Obtain vectors from coordinates
        std::vector<float> d1{
            (atm_2[0] - atm_1[0]), (atm_2[1] - atm_1[1]), (atm_2[2] - atm_1[2])
        };
        std::vector<float> d2{
            (atm_3[0] - atm_2[0]), (atm_3[1] - atm_2[1]), (atm_3[2] - atm_2[2])
        };
        std::vector<float> d3{
            (atm_4[0] - atm_3[0]), (atm_4[1] - atm_3[1]), (atm_4[2] - atm_3[2])
        };
        // 02. Calculate cross product
        std::vector<float> u1 = crossProduct(d1, d2);
        std::vector<float> u2 = crossProduct(d2, d3);

        print3DFloatVec("d1", d1);
        print3DFloatVec("d2", d2);
        print3DFloatVec("d3", d3);
        print3DFloatVec("u1", u1);
        print3DFloatVec("u2", u2);

        // 03. Get cosine theta of u1 & u2
        float cos_torsion = getCosineTheta(u1, u2);
        std::cout << "cos_torsion: " << cos_torsion << std::endl;

        // 04. Calculate arc cosine
        float torsion_angle;
        if (std::isnan(acos(cos_torsion))) {
            if (cos_torsion < 0) {
                torsion_angle = 180.0;
            } else {
                torsion_angle = 0.0;
            }
        } else {
            torsion_angle = acos(cos_torsion) * 180.0 / PI;
        }

        std::cout << "torsion_angle: " << torsion_angle << std::endl;

        // 05. check torsion.xyz.R line 37
        // IMPORTANT: torsion_angle should be minus in specific case;
        // the vector on the plane beta can be represented as
        // cross product of u2 & d2
        std::vector<float> plane_beta_vec = crossProduct(u2, d2);
        std::vector<float> det{
            (u1[0] * plane_beta_vec[0]), (u1[1] * plane_beta_vec[1]),
            (u1[2] * plane_beta_vec[2])
        };

        print3DFloatVec("beta", plane_beta_vec);
        print3DFloatVec("det", det);

        if (det[0] + det[1] + det[2] < 0) {
            torsion_angle = -1 * torsion_angle;
        }
        torsion_vector.push_back(torsion_angle);
        // 06. next index;
        atm_inds[0] += atm_inc;
        atm_inds[1] += atm_inc;
        atm_inds[2] += atm_inc;
        atm_inds[3] += atm_inc;
    }
    return torsion_vector;
}


std::vector<float> getTorsionFromXYZ(
    std::vector<AtomCoordinate> coordinates, int atm_inc
) {
    std::vector< std::vector<float> > coordinate_vector;
    coordinate_vector = extractCoordinates(coordinates);
    std::vector<float> output;
    output = getTorsionFromXYZ(coordinate_vector, atm_inc);
    return output;
}

/**
 * @brief Get the torsion from xyz object
 *
 * @param coordinates a vector of float vectors
 * @param atm_inc an integer
 * @return std::vector<std::vector<float>>
 */
std::vector<float> getTorsionFromXYZ(
    std::vector< std::vector<float> > coordinates, int atm_inc = 1
) {
    std::vector<float> torsion_vector;
    std::vector<int> atm_inds {0,1,2,3};
    while (atm_inds[3] < coordinates.size()) {
        // 00. Get 4 atom coordinates
        std::vector<float> atm_1 = coordinates[atm_inds[0]];
        std::vector<float> atm_2 = coordinates[atm_inds[1]];
        std::vector<float> atm_3 = coordinates[atm_inds[2]];
        std::vector<float> atm_4 = coordinates[atm_inds[3]];
        // 01. Obtain vectors from coordinates
        std::vector<float> d1 {
            (atm_2[0] - atm_1[0]), (atm_2[1] - atm_1[1]), (atm_2[2] - atm_1[2])
        };
        std::vector<float> d2{
            (atm_3[0] - atm_2[0]), (atm_3[1] - atm_2[1]), (atm_3[2] - atm_2[2])
        };
        std::vector<float> d3{
            (atm_4[0] - atm_3[0]), (atm_4[1] - atm_3[1]), (atm_4[2] - atm_3[2])
        };
        // 02. Calculate cross product
        std::vector<float> u1 = crossProduct(d1, d2);
        std::vector<float> u2 = crossProduct(d2, d3);

        // 03. Get cosine theta of u1 & u2
        float cos_torsion = getCosineTheta(u1, u2);
        // 04. Calculate arc cosine
        float torsion_angle;
        if (std::isnan(acos(cos_torsion))) {
            if (cos_torsion < 0) {
                torsion_angle = 180.0;
            } else {
                torsion_angle = 0.0;
            }
        } else {
            torsion_angle = acos(cos_torsion) * 180.0 / PI;
        }

        // 05. check torsion.xyz.R line 37
        // IMPORTANT: torsion_angle should be minus in specific case;
        // the vector on the plane beta can be represented as
        // cross product of u2 & d2
        std::vector<float> plane_beta_vec = crossProduct(u2, d2);
        std::vector<float> det {
            (u1[0] * plane_beta_vec[0]), (u1[1] * plane_beta_vec[1]),
            (u1[2] * plane_beta_vec[2])
        };
        if (det[0] + det[1] + det[2] < 0) {
            torsion_angle = -1 * torsion_angle;
        }
        torsion_vector.push_back(torsion_angle);
        // 06. next index;
        atm_inds[0] += atm_inc;
        atm_inds[1] += atm_inc;
        atm_inds[2] += atm_inc;
        atm_inds[3] += atm_inc;
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
void float3dVectorToDoubleArray(std::vector<float> fv, double output[3]) {
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
void writeTorsionAngles(std::string file_path, std::vector<float> torsion) {
    std::ofstream output;
    output.open(file_path, std::ios_base::out);
    for (float angle : torsion) {
        output << angle << std::endl;
    }
    output.close();
}

std::vector<float> readTorsionAngles(std::string file_path) {
    std::ifstream input_file;
    std::string line;
    std::vector<float> output;
    float angle;
    input_file.open(file_path, std::ios_base::in);
    while(getline(input_file, line)) {
        angle = std::stof(line);
        output.push_back(angle);
    }
    input_file.close();
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
    std::vector<float> torsions, unsigned int n_bits
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
    std::vector<short> encoded_torsions, unsigned int n_bits
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

