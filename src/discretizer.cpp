/**
 * File: discretizer.cpp
 * Project: src
 * Created: 2021-02-05 13:41:54
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     This code is written as part of project "src".
 * ---
 * Last Modified: 2022-06-08 16:53:38
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
#include "discretizer.h"

Discretizer::Discretizer(std::vector<float>& values, unsigned int nb):
    n_bin(nb) {
    // Get min & max
    this->min = *std::min_element(values.begin(), values.end());
    this->max = *std::max_element(values.begin(), values.end());
    // Calculate factors
    this->disc_f = this->n_bin / (this->max - this->min);
    this->cont_f = (this->max - this->min) / this->n_bin;
}

void Discretizer::set_continuous_values(std::vector<float>& values) {
    this->min = *std::min_element(values.begin(), values.end());
    this->max = *std::max_element(values.begin(), values.end());
    // Calculate factors
    this->disc_f = this->n_bin / (this->max - this->min);
    this->cont_f = (this->max - this->min) / this->n_bin;
}

std::vector<unsigned int> Discretizer::discretize(std::vector<float>& continuous_values) {
    std::vector<float>::iterator it;
    unsigned int tmp_disc_value;
    std::vector<unsigned int> discretizedValues;
    for (it = continuous_values.begin(); it != continuous_values.end(); it++) {
        tmp_disc_value = (*it - min) * (this->disc_f);
        discretizedValues.push_back(tmp_disc_value);
    }
    return discretizedValues;
}

unsigned int Discretizer::discretize(float continuous_value) {
    return (continuous_value - this->min) * (this->disc_f);
}

std::vector<float> Discretizer::continuize(std::vector<unsigned int>& discrete_values) {
    std::vector<float> output;
    std::vector<unsigned int>::iterator it;
    float tmp_cont_value;
    for (it = discrete_values.begin(); it != discrete_values.end(); it++) {
        tmp_cont_value = (*it * this->cont_f) + this->min;
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
    fout << "#MIN:" << this->min << std::endl;
    fout << "#MAX:" << this->max << std::endl;
    fout << "#N_BIN:" << this->n_bin << std::endl;
    fout << "#DISC_F:" << this->disc_f << std::endl;
    fout << "#CONT_F:" << this->cont_f << std::endl;
    fout << "ORIGINAL_VALUES,DISCRETIZED_VALUES" << std::endl;
    fout.close();
}

float Discretizer::average_error(std::vector<float>& continuous_values) {
    std::vector<unsigned int> discretized_values = this->discretize(continuous_values);
    std::vector<float> restored = this->continuize(discretized_values);
    float sum = 0;
    for (int i = 0; i < continuous_values.size(); i++) {
        sum += std::abs(continuous_values[i] - restored[i]);
    }
    return sum / continuous_values.size();
}

float Discretizer::max_error(std::vector<float>& continuous_values) {
    std::vector<unsigned int> discretized_values = this->discretize(continuous_values);
    std::vector<float> restored = this->continuize(discretized_values);
    float max = 0;
    for (int i = 0; i < continuous_values.size(); i++) {
        if (std::abs(continuous_values[i] - restored[i]) > max) {
            max = std::abs(continuous_values[i] - restored[i]);
        }
    }
    return max;
}
