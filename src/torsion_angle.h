/**
 * File: torsion_angle.h
 * Project: foldcomp
 * Created: 2021-01-13 10:35:56
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     Functions for torsion angle calculation
 * ---
 * Last Modified: 2022-09-13 15:14:10
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
#pragma once
#include "float3d.h"

#include <cmath>
#include <string> // IWYU pragma: keep
#include <vector>

class AtomCoordinate;

struct TorsionTrace {
    float d2Norm = 0.0f;
    float x = 0.0f;
    float y = 0.0f;
    float angle = 0.0f;
};

static inline TorsionTrace computeTorsionTraceFromCoords(
    const float3d& atm1,
    const float3d& atm2,
    const float3d& atm3,
    const float3d& atm4
) {
    float3d d1{
        atm2.x - atm1.x, atm2.y - atm1.y, atm2.z - atm1.z
    };
    float3d d2{
        atm3.x - atm2.x, atm3.y - atm2.y, atm3.z - atm2.z
    };
    float3d d3{
        atm4.x - atm3.x, atm4.y - atm3.y, atm4.z - atm3.z
    };
    float3d u1{
        std::fmaf(d1.y, d2.z, -(d2.y * d1.z)),
        std::fmaf(d1.z, d2.x, -(d2.z * d1.x)),
        std::fmaf(d1.x, d2.y, -(d2.x * d1.y))
    };
    float3d u2{
        std::fmaf(d2.y, d3.z, -(d3.y * d2.z)),
        std::fmaf(d2.z, d3.x, -(d3.z * d2.x)),
        std::fmaf(d2.x, d3.y, -(d3.x * d2.y))
    };
    float d2Norm = std::sqrt(std::fmaf(d2.x, d2.x, std::fmaf(d2.y, d2.y, d2.z * d2.z)));
    TorsionTrace trace;
    trace.d2Norm = d2Norm;
    if (d2Norm <= FLOAT3D_EPSILON) {
        return trace;
    }
    trace.x = std::fmaf(u1.x, u2.x, std::fmaf(u1.y, u2.y, u1.z * u2.z));
    float3d planeBetaVec{
        std::fmaf(u2.y, d2.z, -(d2.y * u2.z)),
        std::fmaf(u2.z, d2.x, -(d2.z * u2.x)),
        std::fmaf(u2.x, d2.y, -(d2.x * u2.y))
    };
    trace.y = std::fmaf(u1.x, planeBetaVec.x,
                        std::fmaf(u1.y, planeBetaVec.y, u1.z * planeBetaVec.z)) / d2Norm;
    trace.angle = std::atan2(trace.y, trace.x) * 180.0f / static_cast<float>(M_PI);
    return trace;
}

static inline float computeTorsionAngleFromCoords(
    const float3d& atm1,
    const float3d& atm2,
    const float3d& atm3,
    const float3d& atm4
) {
    return computeTorsionTraceFromCoords(atm1, atm2, atm3, atm4).angle;
}

TorsionTrace computeTorsionTraceRuntime(
    const float3d& atm1,
    const float3d& atm2,
    const float3d& atm3,
    const float3d& atm4
);

std::vector<float> getTorsionFromXYZ(
    const std::vector<float3d>& coordinates, int atm_inc
);

std::vector<float> getTorsionFromXYZ(
    const std::vector<AtomCoordinate>& coordinates, int atm_inc
);

std::vector<float> getTorsionFromXYZ(float3d coordinates[4], int atm_inc);

void float3dVectorToDoubleArray(const std::vector<float>& fv, double output[3]);

void writeTorsionAngles(const std::string& file_path, const std::vector<float>& torsion);
std::vector<float> readTorsionAngles(const std::string& file_path);


// 2021-02-02 13:02:02 - encode and decode torsion angles to short

std::vector<short> encodeTorsionAnglesToShort(
    const std::vector<float>& torsions, unsigned int n_bits = 16
);
std::vector<float> decodeEncodedTorsionAngles(
    const std::vector<short>& encoded_torsions, unsigned int n_bits = 16
);
