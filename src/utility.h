/**
 * File: utility.h
 * Project: foldcomp
 * Created: 2021-01-05 14:27:35
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     Utility functions
 * ---
 * Last Modified: 2022-07-20 01:52:32
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
#pragma once
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <fstream>
#include <iostream>
#include <chrono>
#include <sstream>
#include <dirent.h>

const float PI = 3.14159265;

template<typename T>
std::vector<T> vectorSlice(std::vector<T> const& v, int m, int n) {
    auto first = v.cbegin() + m;
    auto last = v.cbegin() + n + 1;

    std::vector<T> vec(first, last);
    return vec;
};

std::vector<float> crossProduct(std::vector<float> v1, std::vector<float> v2);
float norm(std::vector<float> v);
float getCosineTheta(std::vector<float> v1, std::vector<float> v2);

float angle(
    std::vector<float> atm1,
    std::vector<float> atm2,
    std::vector<float> atm3
);

float distance(std::vector<float> atm1, std::vector<float> atm2);



template<typename T>
/**
 * @brief Add the values of map2 to map1
 *
 * @param map1
 * @param map2
 * @return int
 */
int addMap(std::map<std::string, T>* m1, std::map<std::string, T>* m2) {
    size_t m1Size = m1->size();
    if (m1Size == 0) {
        // if m1 is empty, fill in with m2
        for (auto const& x : *m2) {
            (*m1)[x.first] = x.second;
        }
    } else {
        for (auto const& x : *m2) {
            // if the key from map2 is in map1
            if (m1->find(x.first) != m1->end()) {
                // add the values to map1
                (*m1)[x.first] += x.second;
            }
        }
    }
    return 0;
};

//
template <typename T>
int divideMapWithConst(std::map<std::string, T>* m, T c) {
    for (auto const& x: *m) {
        (*m)[x.first] = (x.second / c);
    }
    return 0;
};

template <typename T>
int printMapToFile(std::map<std::string, T>* m, std::string fileName) {
    std::ofstream outFile;
    outFile.open(fileName, std::ios::out);

    for (auto const& x : *m) {
        outFile << x.first << "," << (*m)[x.first] << "\n";
    }
    outFile.close();
    return 0;
};

template <typename T>
void swap(T* a, T* b) {
    T temp = *a;
    *a = *b;
    *b = temp;
};

template <typename T>
void printVector(std::vector<T> v) {
    // Open brackets
    std::cout << "[";
    // Print each element
    for (auto const& x : v) {
        std::cout << x << " ";
    }
    // Close brackets
    std::cout << "]" << std::endl;
};


/**
 * @brief General function to measure the running time of a function
 * https://stackoverflow.com/questions/22387586/measuring-execution-time-of-a-function-in-c
 * @param func
 * @return double
 */
auto static measureRunningTime = [](int n, auto && func, auto&&... params) {
    // get time before function invocation
    const auto& start = std::chrono::high_resolution_clock::now();
    // function invocation using perfect forwarding
    // n is the number of times to run the function; default is 100000
    for (auto i = 0; i < n; ++i) {
        std::forward<decltype(func)>(func)(std::forward<decltype(params)>(params)...);
    }
    // get time after function invocation
    const auto& stop = std::chrono::high_resolution_clock::now();
    return (stop - start) / n;
};

template<class Resolution = std::chrono::milliseconds>
class ExecutionTimer {
public:
    using Clock = std::conditional_t<std::chrono::high_resolution_clock::is_steady,
        std::chrono::high_resolution_clock,
        std::chrono::steady_clock>;
private:
    const Clock::time_point mStart = Clock::now();

public:
    ExecutionTimer() = default;
    ~ExecutionTimer() {
        const auto end = Clock::now();
        std::ostringstream strStream;
        strStream << "Destructor Elapsed: "
            << std::chrono::duration_cast<Resolution>(end - mStart).count()
            << std::endl;
        std::cout << strStream.str() << std::endl;
    }

    inline void stop() {
        const auto end = Clock::now();
        std::ostringstream strStream;
        strStream << "Stop Elapsed: "
            << std::chrono::duration_cast<Resolution>(end - mStart).count()
            << std::endl;
        std::cout << strStream.str() << std::endl;
    }

}; // ExecutionTimer

// Get all files in a directory using dirent.h
std::vector<std::string> getFilesInDirectory(std::string dir);
std::string baseName(std::string const& path);
// One letter code to Three letter code
char getOneLetterCode(std::string threeLetterCode);
std::string getThreeLetterCode(char oneLetterCode);

// Int to One letter code
char convertIntToOneLetterCode(unsigned int aab);
unsigned int convertOneLetterCodeToInt(char oneLetterCode);

// Int to Three letter code
std::string convertIntToThreeLetterCode(unsigned int aab);
unsigned int convertThreeLetterCodeToInt(std::string threeLetterCode);

// CONSTANTS FOR AMINO ACID NAMING
#define AA_ALA_INT 0
#define AA_ALA_STR "ALA"
#define AA_ALA_CHAR 'A'
#define AA_ARG_INT 1
#define AA_ARG_STR "ARG"
#define AA_ARG_CHAR 'R'
#define AA_ASN_INT 2
#define AA_ASN_STR "ASN"
#define AA_ASN_CHAR 'N'
#define AA_ASP_INT 3
#define AA_ASP_STR "ASP"
#define AA_ASP_CHAR 'D'
#define AA_CYS_INT 4
#define AA_CYS_STR "CYS"
#define AA_CYS_CHAR 'C'
#define AA_GLN_INT 5
#define AA_GLN_STR "GLN"
#define AA_GLN_CHAR 'Q'
#define AA_GLU_INT 6
#define AA_GLU_STR "GLU"
#define AA_GLU_CHAR 'E'
#define AA_GLY_INT 7
#define AA_GLY_STR "GLY"
#define AA_GLY_CHAR 'G'
#define AA_HIS_INT 8
#define AA_HIS_STR "HIS"
#define AA_HIS_CHAR 'H'
#define AA_ILE_INT 9
#define AA_ILE_STR "ILE"
#define AA_ILE_CHAR 'I'
#define AA_LEU_INT 10
#define AA_LEU_STR "LEU"
#define AA_LEU_CHAR 'L'
#define AA_LYS_INT 11
#define AA_LYS_STR "LYS"
#define AA_LYS_CHAR 'K'
#define AA_MET_INT 12
#define AA_MET_STR "MET"
#define AA_MET_CHAR 'M'
#define AA_PHE_INT 13
#define AA_PHE_STR "PHE"
#define AA_PHE_CHAR 'F'
#define AA_PRO_INT 14
#define AA_PRO_STR "PRO"
#define AA_PRO_CHAR 'P'
#define AA_SER_INT 15
#define AA_SER_STR "SER"
#define AA_SER_CHAR 'S'
#define AA_THR_INT 16
#define AA_THR_STR "THR"
#define AA_THR_CHAR 'T'
#define AA_TRP_INT 17
#define AA_TRP_STR "TRP"
#define AA_TRP_CHAR 'W'
#define AA_TYR_INT 18
#define AA_TYR_STR "TYR"
#define AA_TYR_CHAR 'Y'
#define AA_VAL_INT 19
#define AA_VAL_STR "VAL"
#define AA_VAL_CHAR 'V'
#define AA_ASX_INT 20
#define AA_ASX_STR "ASX"
#define AA_ASX_CHAR 'B'
#define AA_GLX_INT 21
#define AA_GLX_STR "GLX"
#define AA_GLX_CHAR 'Z'
#define AA_STOP_INT 22
#define AA_STOP_STR "STP"
#define AA_STOP_CHAR '*'
#define AA_UNK_INT 23
#define AA_UNK_STR "UNK"
#define AA_UNK_CHAR 'X'
#define NUM_ALL_AA_CODES 24
#define NUM_VALID_AA_CODES 20