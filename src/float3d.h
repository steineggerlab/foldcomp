#pragma once
#include <cmath>

static constexpr float FLOAT3D_EPSILON = 1e-12f;

struct float3d {
    float3d() : x(0), y(0), z(0) {}
    float3d(float x, float y, float z) : x(x), y(y), z(z) {};
    float x;
    float y;
    float z;
};

/**
 * @brief Return the cross product of two vectors
 *
 * @param v1 A 3d vector of float
 * @param v2 A 3d vector of float
 * @return std::vector<float>
 */
static inline float3d crossProduct(float3d v1, float3d v2) {
    float x = v1.y * v2.z - v2.y * v1.z;
    float y = v1.z * v2.x - v2.z * v1.x;
    float z = v1.x * v2.y - v2.x * v1.y;
    return { x, y, z };
}

static inline float dotProduct(float3d v1, float3d v2) {
    return (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);
}

static inline float normSquared(float3d v) {
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

/**
 * @brief Return the norm of given vector
 *
 * @param v A 3d vector of float
 * @return float
 */
static inline float norm(float3d v) {
    return std::sqrt(normSquared(v));
}

static inline float getCosineTheta(float3d v1, float3d v2) {
    // Calculate inner product of two vectors
    float inner_product = dotProduct(v1, v2);
    float v1_size = normSquared(v1);
    float v2_size = normSquared(v2);
    float denom = std::sqrt(v1_size * v2_size);
    if (denom <= FLOAT3D_EPSILON) {
        return 1.0f;
    }
    float output = inner_product / denom;
    if (output > 1.0f) {
        return 1.0f;
    }
    if (output < -1.0f) {
        return -1.0f;
    }
    return output;
}

static inline float distance(float3d atm1, float3d atm2) {
    float dx = atm1.x - atm2.x;
    float dy = atm1.y - atm2.y;
    float dz = atm1.z - atm2.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

static inline float angle(float3d atm1, float3d atm2, float3d atm3) {
    float3d d1{
        (atm1.x - atm2.x), (atm1.y - atm2.y), (atm1.z - atm2.z)
    };
    float3d d2{
        (atm3.x - atm2.x), (atm3.y - atm2.y), (atm3.z - atm2.z)
    };
    float3d cross = crossProduct(d1, d2);
    float cross_norm = norm(cross);
    float dot = dotProduct(d1, d2);
    float theta = std::atan2(cross_norm, dot) * 180.0f / static_cast<float>(M_PI);
    return theta;
}
