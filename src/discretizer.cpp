/**
 * File: discretizer.cpp
 * Project: foldcomp
 * Created: 2021-02-05 13:41:54
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     Functions for discretizing float values and restoring them
 * ---
 * Last Modified: 2022-12-09 15:42:34
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
#include "discretizer.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream> // IWYU pragma: keep
#include <iostream>

namespace {

void setDiscretizerRange(Discretizer& discretizer, float minValue, float maxValue) {
    discretizer.min = minValue;
    discretizer.max = maxValue;
    const float range = discretizer.max - discretizer.min;
    if (discretizer.n_bin == 0 || std::fabs(range) <= 1e-12f) {
        discretizer.disc_f = 0.0f;
        discretizer.cont_f = 0.0f;
        return;
    }
    discretizer.disc_f = static_cast<float>(discretizer.n_bin) / range;
    discretizer.cont_f = range / static_cast<float>(discretizer.n_bin);
}

unsigned int discretizeValue(const Discretizer& discretizer, float continuousValue) {
    if (discretizer.n_bin == 0 || discretizer.disc_f == 0.0f) {
        return 0;
    }
    float scaled = (continuousValue - discretizer.min) * discretizer.disc_f;
    long rounded = std::lround(scaled);
    if (rounded < 0) {
        return 0;
    }
    if (rounded > static_cast<long>(discretizer.n_bin)) {
        return discretizer.n_bin;
    }
    return static_cast<unsigned int>(rounded);
}

}

Discretizer::Discretizer(const std::vector<float>& values, unsigned int nb):
    n_bin(nb) {
    if (values.size() == 0) {
        min = 0.0f;
        max = 0.0f;
        disc_f = 0.0f;
        cont_f = 0.0f;
        return;
    }
    setDiscretizerRange(*this,
                        *std::min_element(values.begin(), values.end()),
                        *std::max_element(values.begin(), values.end()));
}

void Discretizer::set_continuous_values(const std::vector<float>& values) {
    if (values.empty()) {
        setDiscretizerRange(*this, 0.0f, 0.0f);
        return;
    }
    setDiscretizerRange(*this,
                        *std::min_element(values.begin(), values.end()),
                        *std::max_element(values.begin(), values.end()));
}

std::vector<unsigned int> Discretizer::discretize(const std::vector<float>& continuous_values) {
    std::vector<unsigned int> discretizedValues;
    discretizedValues.reserve(continuous_values.size());
    for (float value : continuous_values) {
        discretizedValues.push_back(discretizeValue(*this, value));
    }
    return discretizedValues;
}

unsigned int Discretizer::discretize(float continuous_value) {
    return discretizeValue(*this, continuous_value);
}

std::vector<float> Discretizer::continuize(const std::vector<unsigned int>& discrete_values) {
    std::vector<unsigned int>::const_iterator it;
    std::vector<float> output;
    output.reserve(discrete_values.size());
    for (it = discrete_values.cbegin(); it != discrete_values.cend(); it++) {
        float tmp_cont_value = (*it * this->cont_f) + this->min;
        output.push_back(tmp_cont_value);
    }
    return output;
}

float Discretizer::continuize(unsigned int discrete_value) {
    return (discrete_value * this->cont_f) + this->min;
}

DiscParams Discretizer::get_param() {
    DiscParams params;
    params.min = this->min;
    params.max = this->max;
    params.n_bin = this->n_bin;
    params.disc_f = this->disc_f;
    params.cont_f = this->cont_f;
    return params;
}

// Methods for tests

void Discretizer::print() {
    std::cout << "MIN: " << this->min << std::endl;
    std::cout << "MAX: " << this->max << std::endl;
    std::cout << "N_BIN: " << this->n_bin << std::endl;
    std::cout << "DISC_F: " << this->disc_f << std::endl;
    std::cout << "CONT_F: " << this->cont_f << std::endl;
}

void Discretizer::write_to_file(std::string filename) {
    std::ofstream fout(filename);
    fout << "#MIN:" << this->min << "\n";
    fout << "#MAX:" << this->max << "\n";
    fout << "#N_BIN:" << this->n_bin << "\n";
    fout << "#DISC_F:" << this->disc_f << "\n";
    fout << "#CONT_F:" << this->cont_f << "\n";
    fout << "ORIGINAL_VALUES,DISCRETIZED_VALUES\n";
    fout.close();
}

float Discretizer::average_error(const std::vector<float>& continuous_values) {
    std::vector<unsigned int> discretized_values = this->discretize(continuous_values);
    std::vector<float> restored = this->continuize(discretized_values);
    float sum = 0;
    for (size_t i = 0; i < continuous_values.size(); i++) {
        sum += std::abs(continuous_values[i] - restored[i]);
    }
    return sum / continuous_values.size();
}

float Discretizer::max_error(const std::vector<float>& continuous_values) {
    std::vector<unsigned int> discretized_values = this->discretize(continuous_values);
    std::vector<float> restored = this->continuize(discretized_values);
    float max = 0;
    for (size_t i = 0; i < continuous_values.size(); i++) {
        if (std::abs(continuous_values[i] - restored[i]) > max) {
            max = std::abs(continuous_values[i] - restored[i]);
        }
    }
    return max;
}
