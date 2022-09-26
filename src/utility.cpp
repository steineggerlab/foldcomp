/**
 * File: utility.cpp
 * Project: foldcomp
 * Created: 2021-01-05 14:29:04
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     Utility functions
 * ---
 * Last Modified: 2022-09-20 11:51:28
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */

// TODO: TRIM WHITE SPACES FROM PDB LINE

#include "utility.h"

/**
 * @brief Return the cross product of two vectors
 *
 * @param v1 A 3d vector of float
 * @param v2 A 3d vector of float
 * @return std::vector<float>
 */
std::vector<float> crossProduct(std::vector<float> v1, std::vector<float> v2) {
    // TODO: Check the length of v1, v2 to be 3
    float x = (v1[1] * v2[2]) - (v2[1] * v1[2]);
    float y = (v1[2] * v2[0]) - (v2[2] * v1[0]);
    float z = (v1[0] * v2[1]) - (v2[0] * v1[1]);
    // Calculate
    std::vector<float> product{ x, y, z };
    return product;
}

/**
 * @brief Return the norm of given vector
 *
 * @param v A 3d vector of float
 * @return float
 */
float norm(std::vector<float> v) {
    return sqrt(pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2));
}

float getCosineTheta(std::vector<float> v1, std::vector<float> v2) {
    float output;
    // TODO: Check the length of v1, v2 to be 3
    // Calculate inner product of two vectors
    float inner_product = (v1[0] * v2[0]) + (v1[1] * v2[1]) + (v1[2] * v2[2]);
    float v1_size = pow(v1[0], 2) + pow(v1[1], 2) + pow(v1[2], 2);
    float v2_size = pow(v2[0], 2) + pow(v2[1], 2) + pow(v2[2], 2);
    output = inner_product / sqrt(v1_size * v2_size);
    return output;
}

float distance(std::vector<float> atm1, std::vector<float> atm2) {
    float output = 0.0;
    output = sqrt(
        (pow(atm1[0] - atm2[0], 2) +
            pow(atm1[1] - atm2[1], 2) +
            pow(atm1[2] - atm2[2], 2))
    );
    return output;
}

float angle(
    std::vector<float> atm1,
    std::vector<float> atm2,
    std::vector<float> atm3
) {
    std::vector<float> d1{
        (atm1[0] - atm2[0]), (atm1[1] - atm2[1]), (atm1[2] - atm2[2])
    };
    std::vector<float> d2{
        (atm3[0] - atm2[0]), (atm3[1] - atm2[1]), (atm3[2] - atm2[2])
    };
    float cos_theta = getCosineTheta(d1, d2);
    float theta = acos(cos_theta) * 180.0 / PI;
    return theta;
}

#ifdef FOLDCOMP_EXECUTABLE
#include <dirent.h>
// Get all files in a directory using dirent.h
std::vector<std::string> getFilesInDirectory(std::string dir) {
    std::vector<std::string> files;
    DIR* dp;
    struct dirent* dirp;
    if ((dp = opendir(dir.c_str())) == nullptr) {
        std::cout << "Error(" << errno << ") opening " << dir << std::endl;
        return files;
    }

    while ((dirp = readdir(dp)) != nullptr) {
        // Skip directories
        if (dirp->d_type == DT_DIR) {
            continue;
        }
        files.push_back(std::string(dirp->d_name));
    }
    closedir(dp);
    return files;
};
#endif

std::string baseName(std::string const path) {
    return path.substr(path.find_last_of("/\\") + 1);
}

std::string getFileWithoutExt(std::string& file) {
    size_t extStart = file.find_last_of('.');
    return extStart == std::string::npos ? file : file.substr(0, extStart);
}

bool stringEndsWith(const std::string& suffix, const std::string& str) {
    if (str.length() < suffix.length()) {
        return false;
    }

    return (!str.compare(str.length() - suffix.length(), suffix.length(), suffix));
}

bool stringStartsWith(const std::string& prefix, const std::string& str, const size_t offset) {
    if (str.length() < prefix.length()) {
        return false;
    }
    return (!str.compare(offset, prefix.length(), prefix));
}

std::vector<std::string> stringSplit(const std::string& str, const std::string& sep) {
    std::vector<std::string> arr;

    char* cstr = strdup(str.c_str());
    const char* csep = sep.c_str();
    char* rest;
    char* current = strtok_r(cstr, csep, &rest);
    while (current != NULL) {
        arr.emplace_back(current);
        current = strtok_r(NULL, csep, &rest);
    }
    free(cstr);

    return arr;
}

/*  FUNCTIONS FOR HANDLING AMINO ACID NAME  */

char getOneLetterCode(std::string threeLetterCode) {
    char aa;
    if (threeLetterCode == AA_ALA_STR) {
        aa = AA_ALA_CHAR;
    } else if (threeLetterCode == AA_ARG_STR) {
        aa = AA_ARG_CHAR;
    } else if (threeLetterCode == AA_ASN_STR) {
        aa = AA_ASN_CHAR;
    } else if (threeLetterCode == AA_ASP_STR) {
        aa = AA_ASP_CHAR;
    } else if (threeLetterCode == AA_CYS_STR) {
        aa = AA_CYS_CHAR;
    } else if (threeLetterCode == AA_GLN_STR) {
        aa = AA_GLN_CHAR;
    } else if (threeLetterCode == AA_GLU_STR) {
        aa = AA_GLU_CHAR;
    } else if (threeLetterCode == AA_GLY_STR) {
        aa = AA_GLY_CHAR;
    } else if (threeLetterCode == AA_HIS_STR) {
        aa = AA_HIS_CHAR;
    } else if (threeLetterCode == AA_ILE_STR) {
        aa = AA_ILE_CHAR;
    } else if (threeLetterCode == AA_LEU_STR) {
        aa = AA_LEU_CHAR;
    } else if (threeLetterCode == AA_LYS_STR) {
        aa = AA_LYS_CHAR;
    } else if (threeLetterCode == AA_MET_STR) {
        aa = AA_MET_CHAR;
    } else if (threeLetterCode == AA_PHE_STR) {
        aa = AA_PHE_CHAR;
    } else if (threeLetterCode == AA_PRO_STR) {
        aa = AA_PRO_CHAR;
    } else if (threeLetterCode == AA_SER_STR) {
        aa = AA_SER_CHAR;
    } else if (threeLetterCode == AA_THR_STR) {
        aa = AA_THR_CHAR;
    } else if (threeLetterCode == AA_TRP_STR) {
        aa = AA_TRP_CHAR;
    } else if (threeLetterCode == AA_TYR_STR) {
        aa = AA_TYR_CHAR;
    } else if (threeLetterCode == AA_VAL_STR) {
        aa = AA_VAL_CHAR;
    } else if (threeLetterCode == AA_ASX_STR) {
        aa = AA_ASX_CHAR;
    } else if (threeLetterCode == AA_GLX_STR) {
        aa = AA_GLX_CHAR;
    } else if (threeLetterCode == AA_STOP_STR) {
        aa = AA_STOP_CHAR;
    } else if (threeLetterCode == AA_UNK_STR) {
        aa = AA_UNK_CHAR;
    } else {
        aa = AA_UNK_CHAR;
    }
    return aa;
}

std::string getThreeLetterCode(char oneLetterCode) {
    std::string threeLetterCode;
    if (oneLetterCode == AA_ALA_CHAR) {
        threeLetterCode = AA_ALA_STR;
    } else if (oneLetterCode == AA_ARG_CHAR) {
        threeLetterCode = AA_ARG_STR;
    } else if (oneLetterCode == AA_ASN_CHAR) {
        threeLetterCode = AA_ASN_STR;
    } else if (oneLetterCode == AA_ASP_CHAR) {
        threeLetterCode = AA_ASP_STR;
    } else if (oneLetterCode == AA_CYS_CHAR) {
        threeLetterCode = AA_CYS_STR;
    } else if (oneLetterCode == AA_GLN_CHAR) {
        threeLetterCode = AA_GLN_STR;
    } else if (oneLetterCode == AA_GLU_CHAR) {
        threeLetterCode = AA_GLU_STR;
    } else if (oneLetterCode == AA_GLY_CHAR) {
        threeLetterCode = AA_GLY_STR;
    } else if (oneLetterCode == AA_HIS_CHAR) {
        threeLetterCode = AA_HIS_STR;
    } else if (oneLetterCode == AA_ILE_CHAR) {
        threeLetterCode = AA_ILE_STR;
    } else if (oneLetterCode == AA_LEU_CHAR) {
        threeLetterCode = AA_LEU_STR;
    } else if (oneLetterCode == AA_LYS_CHAR) {
        threeLetterCode = AA_LYS_STR;
    } else if (oneLetterCode == AA_MET_CHAR) {
        threeLetterCode = AA_MET_STR;
    } else if (oneLetterCode == AA_PHE_CHAR) {
        threeLetterCode = AA_PHE_STR;
    } else if (oneLetterCode == AA_PRO_CHAR) {
        threeLetterCode = AA_PRO_STR;
    } else if (oneLetterCode == AA_SER_CHAR) {
        threeLetterCode = AA_SER_STR;
    } else if   (oneLetterCode == AA_THR_CHAR) {
        threeLetterCode = AA_THR_STR;
    } else if (oneLetterCode == AA_TRP_CHAR) {
        threeLetterCode = AA_TRP_STR;
    } else if (oneLetterCode == AA_TYR_CHAR) {
        threeLetterCode = AA_TYR_STR;
    } else if (oneLetterCode == AA_VAL_CHAR) {
        threeLetterCode = AA_VAL_STR;
    } else if (oneLetterCode == AA_ASX_CHAR) {
        threeLetterCode = AA_ASX_STR;
    } else if (oneLetterCode == AA_GLX_CHAR) {
        threeLetterCode = AA_GLX_STR;
    } else if (oneLetterCode == AA_STOP_CHAR) {
        threeLetterCode = AA_STOP_STR;
    } else if (oneLetterCode == AA_UNK_CHAR) {
        threeLetterCode = AA_UNK_STR;
    } else {
        threeLetterCode = AA_UNK_STR;
    }
    return threeLetterCode;
}

/**
 * @brief Read 5-bit encoded residue (unsigned integer) and return
 *       corresponding amino acid (char)
 *
 * @param aab
 * @return char
 */
char convertIntToOneLetterCode(unsigned int aab) {
    char aa;
    switch (aab) {
        case AA_ALA_INT:
            aa = AA_ALA_CHAR;
            break;
        case AA_ARG_INT:
            aa = AA_ARG_CHAR;
            break;
        case AA_ASN_INT:
            aa = AA_ASN_CHAR;
            break;
        case AA_ASP_INT:
            aa = AA_ASP_CHAR;
            break;
        case AA_CYS_INT:
            aa = AA_CYS_CHAR;
            break;
        case AA_GLN_INT:
            aa = AA_GLN_CHAR;
            break;
        case AA_GLU_INT:
            aa = AA_GLU_CHAR;
            break;
        case AA_GLY_INT:
            aa = AA_GLY_CHAR;
            break;
        case AA_HIS_INT:
            aa = AA_HIS_CHAR;
            break;
        case AA_ILE_INT:
            aa = AA_ILE_CHAR;
            break;
        case AA_LEU_INT:
            aa = AA_LEU_CHAR;
            break;
        case AA_LYS_INT:
            aa = AA_LYS_CHAR;
            break;
        case AA_MET_INT:
            aa = AA_MET_CHAR;
            break;
        case AA_PHE_INT:
            aa = AA_PHE_CHAR;
            break;
        case AA_PRO_INT:
            aa = AA_PRO_CHAR;
            break;
        case AA_SER_INT:
            aa = AA_SER_CHAR;
            break;
        case AA_THR_INT:
            aa = AA_THR_CHAR;
            break;
        case AA_TRP_INT:
            aa = AA_TRP_CHAR;
            break;
        case AA_TYR_INT:
            aa = AA_TYR_CHAR;
            break;
        case AA_VAL_INT:
            aa = AA_VAL_CHAR;
            break;
        case AA_ASX_INT:
            aa = AA_ASX_CHAR;
            break;
        case AA_GLX_INT:
            aa = AA_GLX_CHAR;
            break;
        case AA_STOP_INT:
            aa = AA_STOP_CHAR;
            break;
        case AA_UNK_INT:
            aa = AA_UNK_CHAR;
            break;
        default:
            aa = AA_UNK_CHAR;
            break;
    }
    return aa;
}

unsigned int convertOneLetterCodeToInt(char oneLetterCode) {
    unsigned int output;
    switch (oneLetterCode) {
        case AA_ALA_CHAR:
            output = AA_ALA_INT;
            break;
        case AA_ARG_CHAR:
            output = AA_ARG_INT;
            break;
        case AA_ASN_CHAR:
            output = AA_ASN_INT;
            break;
        case AA_ASP_CHAR:
            output = AA_ASP_INT;
            break;
        case AA_CYS_CHAR:
            output = AA_CYS_INT;
            break;
        case AA_GLN_CHAR:
            output = AA_GLN_INT;
            break;
        case AA_GLU_CHAR:
            output = AA_GLU_INT;
            break;
        case AA_GLY_CHAR:
            output = AA_GLY_INT;
            break;
        case AA_HIS_CHAR:
            output = AA_HIS_INT;
            break;
        case AA_ILE_CHAR:
            output = AA_ILE_INT;
            break;
        case AA_LEU_CHAR:
            output = AA_LEU_INT;
            break;
        case AA_LYS_CHAR:
            output = AA_LYS_INT;
            break;
        case AA_MET_CHAR:
            output = AA_MET_INT;
            break;
        case AA_PHE_CHAR:
            output = AA_PHE_INT;
            break;
        case AA_PRO_CHAR:
            output = AA_PRO_INT;
            break;
        case AA_SER_CHAR:
            output = AA_SER_INT;
            break;
        case AA_THR_CHAR:
            output = AA_THR_INT;
            break;
        case AA_TRP_CHAR:
            output = AA_TRP_INT;
            break;
        case AA_TYR_CHAR:
            output = AA_TYR_INT;
            break;
        case AA_VAL_CHAR:
            output = AA_VAL_INT;
            break;
        case AA_ASX_CHAR:
            output = AA_ASX_INT;
            break;
        case AA_GLX_CHAR:
            output = AA_GLX_INT;
            break;
        case AA_STOP_CHAR:
            output = AA_STOP_INT;
            break;
        case AA_UNK_CHAR:
            output = AA_UNK_INT;
            break;
        default:
            output = AA_UNK_INT;
            break;
    }
    return output;
}

std::string convertIntToThreeLetterCode(unsigned int aab) {
    std::string threeLetterCode;
    switch (aab) {
        case AA_ALA_INT:
            threeLetterCode = AA_ALA_STR;
            break;
        case AA_ARG_INT:
            threeLetterCode = AA_ARG_STR;
            break;
        case AA_ASN_INT:
            threeLetterCode = AA_ASN_STR;
            break;
        case AA_ASP_INT:
            threeLetterCode = AA_ASP_STR;
            break;
        case AA_CYS_INT:
            threeLetterCode = AA_CYS_STR;
            break;
        case AA_GLN_INT:
            threeLetterCode = AA_GLN_STR;
            break;
        case AA_GLU_INT:
            threeLetterCode = AA_GLU_STR;
            break;
        case AA_GLY_INT:
            threeLetterCode = AA_GLY_STR;
            break;
        case AA_HIS_INT:
            threeLetterCode = AA_HIS_STR;
            break;
        case AA_ILE_INT:
            threeLetterCode = AA_ILE_STR;
            break;
        case AA_LEU_INT:
            threeLetterCode = AA_LEU_STR;
            break;
        case AA_LYS_INT:
            threeLetterCode = AA_LYS_STR;
            break;
        case AA_MET_INT:
            threeLetterCode = AA_MET_STR;
            break;
        case AA_PHE_INT:
            threeLetterCode = AA_PHE_STR;
            break;
        case AA_PRO_INT:
            threeLetterCode = AA_PRO_STR;
            break;
        case AA_SER_INT:
            threeLetterCode = AA_SER_STR;
            break;
        case AA_THR_INT:
            threeLetterCode = AA_THR_STR;
            break;
        case AA_TRP_INT:
            threeLetterCode = AA_TRP_STR;
            break;
        case AA_TYR_INT:
            threeLetterCode = AA_TYR_STR;
            break;
        case AA_VAL_INT:
            threeLetterCode = AA_VAL_STR;
            break;
        case AA_ASX_INT:
            threeLetterCode = AA_ASX_STR;
            break;
        case AA_GLX_INT:
            threeLetterCode = AA_GLX_STR;
            break;
        case AA_STOP_INT:
            threeLetterCode = AA_STOP_STR;
            break;
        case AA_UNK_INT:
            threeLetterCode = AA_UNK_STR;
            break;
        default:
            threeLetterCode = AA_UNK_STR;
            break;
    }
    return threeLetterCode;
}

unsigned int convertThreeLetterCodeToInt(std::string threeLetterCode) {
    unsigned int output;
    if (threeLetterCode == AA_ALA_STR) {
        output = AA_ALA_INT;
    } else if (threeLetterCode == AA_ARG_STR) {
        output = AA_ARG_INT;
    } else if (threeLetterCode == AA_ASN_STR) {
        output = AA_ASN_INT;
    } else if (threeLetterCode == AA_ASP_STR) {
        output = AA_ASP_INT;
    } else if (threeLetterCode == AA_CYS_STR) {
        output = AA_CYS_INT;
    } else if (threeLetterCode == AA_GLN_STR) {
        output = AA_GLN_INT;
    } else if (threeLetterCode == AA_GLU_STR) {
        output = AA_GLU_INT;
    } else if (threeLetterCode == AA_GLY_STR) {
        output = AA_GLY_INT;
    } else if (threeLetterCode == AA_HIS_STR) {
        output = AA_HIS_INT;
    } else if (threeLetterCode == AA_ILE_STR) {
        output = AA_ILE_INT;
    } else if (threeLetterCode == AA_LEU_STR) {
        output = AA_LEU_INT;
    } else if (threeLetterCode == AA_LYS_STR) {
        output = AA_LYS_INT;
    } else if (threeLetterCode == AA_MET_STR) {
        output = AA_MET_INT;
    } else if (threeLetterCode == AA_PHE_STR) {
        output = AA_PHE_INT;
    } else if (threeLetterCode == AA_PRO_STR) {
        output = AA_PRO_INT;
    } else if (threeLetterCode == AA_SER_STR) {
        output = AA_SER_INT;
    } else if (threeLetterCode == AA_THR_STR) {
        output = AA_THR_INT;
    } else if (threeLetterCode == AA_TRP_STR) {
        output = AA_TRP_INT;
    } else if (threeLetterCode == AA_TYR_STR) {
        output = AA_TYR_INT;
    } else if (threeLetterCode == AA_VAL_STR) {
        output = AA_VAL_INT;
    } else if (threeLetterCode == AA_ASX_STR) {
        output = AA_ASX_INT;
    } else if (threeLetterCode == AA_GLX_STR) {
        output = AA_GLX_INT;
    } else if (threeLetterCode == AA_STOP_STR) {
        output = AA_STOP_INT;
    } else if (threeLetterCode == AA_UNK_STR) {
        output = AA_UNK_INT;
    } else {
        output = AA_UNK_INT;
    }
    return output;
}
