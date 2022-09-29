#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: rmsd_from_pdbs.py
# Project: foldcomp
# Description:
#     A script to calculate the RMSD between two PDB files.
# Usage:
#     python rmsd_from_pdbs.py <pdb_path1> <pdb_path2>
# ------------------------------------------------------------
# Created: 2022-07-04 12:20:45
# Author: Hyunbin Kim (khb7840@gmail.com)
# ------------------------------------------------------------
# Last Modified: 2022-07-04 12:46:56
# Modified By: Hyunbin Kim (khb7840@gmail.com)
# ------------------------------------------------------------
# Copyright © 2022 Hyunbin Kim, All rights reserved

from Bio import PDB
from Bio.PDB import QCPSuperimposer
import numpy as np
import sys


def calculate_rmsd(pdb_path1, pdb_path2, mode="CA"):

    with open(pdb_path1, "r") as f:
        with open(pdb_path2, "r") as f2:
            atom1_list = []
            atom2_list = []
            for line in f:
                if line[0:4] == "ATOM":
                    if mode == "CA":
                        if line[13:15] == "CA":
                            atom1_list.append(
                                np.array(
                                    [
                                        float(line[30:38]),
                                        float(line[38:46]),
                                        float(line[46:54]),
                                    ]
                                )
                            )
                    elif mode == "Backbone":
                        if line[13:15] in ["N", "CA", "C"]:
                            atom1_list.append(
                                np.array(
                                    [
                                        float(line[30:38]),
                                        float(line[38:46]),
                                        float(line[46:54]),
                                    ]
                                )
                            )
                    elif mode == "All":
                        atom1_list.append(
                            np.array(
                                [
                                    float(line[30:38]),
                                    float(line[38:46]),
                                    float(line[46:54]),
                                ]
                            )
                        )
            for line in f2:
                if line[0:4] == "ATOM":
                    if mode == "CA":
                        if line[13:15] == "CA":
                            atom2_list.append(
                                np.array(
                                    [
                                        float(line[30:38]),
                                        float(line[38:46]),
                                        float(line[46:54]),
                                    ]
                                )
                            )
                    elif mode == "Backbone":
                        if line[13:15] in ["N", "CA", "C"]:
                            atom2_list.append(
                                np.array(
                                    [
                                        float(line[30:38]),
                                        float(line[38:46]),
                                        float(line[46:54]),
                                    ]
                                )
                            )
                    elif mode == "All":
                        atom2_list.append(
                            np.array(
                                [
                                    float(line[30:38]),
                                    float(line[38:46]),
                                    float(line[46:54]),
                                ]
                            )
                        )

    sup = QCPSuperimposer.QCPSuperimposer()

    atom1_list = np.array(atom1_list)
    atom2_list = np.array(atom2_list)

    sup.set(atom1_list, atom2_list)
    sup.run()
    rmsd = sup.get_rms()
    return rmsd


def main():
    path1 = sys.argv[1]
    path2 = sys.argv[2]
    backbone_rmsd = calculate_rmsd(path1, path2, "Backbone")
    all_rmsd = calculate_rmsd(path1, path2, "All")
    ca_rmsd = calculate_rmsd(path1, path2, "CA")
    print("C-alpha RMSD: %.3f" % ca_rmsd)
    print("Backbone RMSD: %.3f" % backbone_rmsd)
    print("All RMSD: %.3f" % all_rmsd)


if __name__ == "__main__":
    main()
