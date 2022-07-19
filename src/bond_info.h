/**
 * File: bond_info.h
 * Project: foldcomp
 * Created: 2021-01-26 14:37:08
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     This file contains informations about bonds in protein structures.
 * ---
 * Last Modified: 2022-07-20 01:58:44
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */

#pragma once

#include <string>
#include <map>

class BondInfo {
private:
    /* data */
public:
    BondInfo(/* args */);
    ~BondInfo();

    std::map<std::string, float> aminoAcidBondLengths() {
        std::map<std::string, float> output {
            {"N_TO_CA", 1.46}, {"CA_TO_C", 1.52}, {"C_TO_N", 1.33}
        };
        return output;
    }

    /**
     * @brief Amino-acid specific bond angles.
     * (Pre-defined constants from PeptideBuilder)
     *
     * @return std::map<std::string, std::map<std::string, float> >
     */
    std::map<std::string, std::map<std::string, float> > brafAminoAcidBondAngles(){
        std::map< std::string, std::map<std::string, float> > output;

        // This values are copied from PeptideBuilder
        // C_TO_N & N_TO_CA seems to be constant
        output["ALA"] = {{"CA_TO_C", 111.14246875},
                         {"C_TO_N", 117.21509375},
                         {"N_TO_CA", 121.0790625}};
        output["ARG"] = {{"CA_TO_C", 111.14375},
                         {"C_TO_N", 117.259928571429},
                         {"N_TO_CA", 121.589678571429}};
        output["ASN"] = {{"CA_TO_C", 111.3495},
                         {"C_TO_N", 117.65195},
                         {"N_TO_CA", 121.58975}};
        output["ASP"] = {{"CA_TO_C", 111.193161290323},
                         {"C_TO_N", 116.918935483871},
                         {"N_TO_CA", 121.293193548387}};
        output["CYS"] = {{"CA_TO_C", 110.09825},
                         {"C_TO_N", 116.874875},
                         {"N_TO_CA", 121.371875}};
        output["GLN"] = {{"CA_TO_C", 111.463533333333},
                         {"C_TO_N", 117.168931034483},
                         {"N_TO_CA", 121.0669}};
        output["GLU"] = {{"CA_TO_C", 110.810513513514},
                         {"C_TO_N", 116.901631578947},
                         {"N_TO_CA", 121.225837837838}};
        output["GLY"] = {{"CA_TO_C", 113.057534883721},
                         {"C_TO_N", 116.560558139535},
                         {"N_TO_CA", 121.06223255814}};
        output["HIS"] = {{"CA_TO_C", 110.647647058824},
                         {"C_TO_N", 116.874529411765},
                         {"N_TO_CA", 121.465235294118}};
        output["ILE"] = {{"CA_TO_C", 109.595130434783},
                         {"C_TO_N", 117.084260869565},
                         {"N_TO_CA", 121.780717391304}};
        output["LEU"] = {{"CA_TO_C", 111.10562295082},
                         {"C_TO_N", 117.091262295082},
                         {"N_TO_CA", 121.130540983607}};
        output["LYS"] = {{"CA_TO_C", 110.861119047619},
                         {"C_TO_N", 117.27519047619},
                         {"N_TO_CA", 121.480166666667}};
        output["MET"] = {{"CA_TO_C", 111.076},
                         {"C_TO_N", 113.171995238095},
                         {"N_TO_CA", 121.297476190476}};
        output["PHE"] = {{"CA_TO_C", 110.71925},
                         {"C_TO_N", 116.78825},
                         {"N_TO_CA", 121.350541666667}};
        output["PRO"] = {{"CA_TO_C", 112.703357142857},
                         {"C_TO_N", 116.532678571429},
                         {"N_TO_CA", 121.274642857143}};
        output["SER"] = {{"CA_TO_C", 110.630595238095},
                         {"C_TO_N", 116.829119047619},
                         {"N_TO_CA", 121.735073170732}};
        output["THR"] = {{"CA_TO_C", 110.790294117647},
                         {"C_TO_N", 117.127705882353},
                         {"N_TO_CA", 121.702647058824}};
        output["TRP"] = {{"CA_TO_C", 110.835714285714},
                         {"C_TO_N", 116.895571428571},
                         {"N_TO_CA", 121.944714285714}};
        output["TYR"] = {{"CA_TO_C", 111.080058823529},
                         {"C_TO_N", 117.017117647059},
                         {"N_TO_CA", 121.753058823529}};
        output["VAL"] = {{"CA_TO_C", 109.875026315789},
                         {"C_TO_N", 116.857947368421},
                         {"N_TO_CA", 121.607815789474}};
        return output;
    }


    /**
     * @brief Amino-acid specific bond angles.
     * (Pre-defined constants from PeptideBuilder)
     *
     * @return std::map<std::string, std::map<std::string, float> >
     */
    std::map<std::string, std::map<std::string, float> > aminoAcidBondAngles() {
      std::map<std::string, std::map<std::string, float> > output;

      // This values are copied from PeptideBuilder
      // C_TO_N & N_TO_CA seems to be constant
      output["GLY"] = {
          {"CA_TO_C", 116.5605}, {"C_TO_N", 121.0622}, {"N_TO_CA", 113.0575}};
      output["ALA"] = {
          {"CA_TO_C", 111.068}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["SER"] = {
          {"CA_TO_C", 111.2812}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["CYS"] = {
          {"CA_TO_C", 110.8856}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["VAL"] = {
          {"CA_TO_C", 111.068}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["ILE"] = {
          {"CA_TO_C", 109.7202}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["LEU"] = {
          {"CA_TO_C", 110.8652}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["THR"] = {
          {"CA_TO_C", 110.7014}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["ARG"] = {
          {"CA_TO_C", 110.98}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["LYS"] = {
          {"CA_TO_C", 111.08}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["ASP"] = {
          {"CA_TO_C", 111.03}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["ASN"] = {
          {"CA_TO_C", 111.5}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["GLU"] = {
          {"CA_TO_C", 111.1703}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["GLN"] = {
          {"CA_TO_C", 111.0849}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["MET"] = {
          {"CA_TO_C", 110.9416}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["HIS"] = {
          {"CA_TO_C", 111.0859}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["PRO"] = {
          {"CA_TO_C", 112.7499}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["PHE"] = {
          {"CA_TO_C", 110.7528}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["TYR"] = {
          {"CA_TO_C", 110.9288}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["TRP"] = {
          {"CA_TO_C", 110.8914}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};

      return output;
    }

    /**
     * @brief Amino-acid specific bond angles.
     * (Pre-defined constants from PeptideBuilder)
     *
     * @return std::map<std::string, std::map<std::string, float> >
     */
    std::map<std::string, std::map<std::string, float> > brafAminoAcidBondLengths() {
      std::map<std::string, std::map<std::string, float> > output;

      // This values are copied from PeptideBuilder
      // C_TO_N & N_TO_CA seems to be constant
      output["GLY"] = {
          {"CA_TO_C", 110.8914}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["ALA"] = {
          {"CA_TO_C", 111.068}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["SER"] = {
          {"CA_TO_C", 111.2812}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["CYS"] = {
          {"CA_TO_C", 110.8856}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["VAL"] = {
          {"CA_TO_C", 111.068}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["ILE"] = {
          {"CA_TO_C", 109.7202}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["LEU"] = {
          {"CA_TO_C", 110.8652}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["THR"] = {
          {"CA_TO_C", 110.7014}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["ARG"] = {
          {"CA_TO_C", 110.98}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["LYS"] = {
          {"CA_TO_C", 111.08}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["ASP"] = {
          {"CA_TO_C", 111.03}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["ASN"] = {
          {"CA_TO_C", 111.5}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["GLU"] = {
          {"CA_TO_C", 111.1703}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["GLN"] = {
          {"CA_TO_C", 111.0849}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["MET"] = {
          {"CA_TO_C", 110.9416}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["HIS"] = {
          {"CA_TO_C", 111.0859}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["PRO"] = {
          {"CA_TO_C", 112.7499}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["PHE"] = {
          {"CA_TO_C", 110.7528}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["TYR"] = {
          {"CA_TO_C", 110.9288}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};
      output["TRP"] = {
          {"CA_TO_C", 110.8914}, {"C_TO_N", 116.643}, {"N_TO_CA", 121.3822}};

      return output;
    }
};

BondInfo::BondInfo(/* args */) {}

BondInfo::~BondInfo() {}
