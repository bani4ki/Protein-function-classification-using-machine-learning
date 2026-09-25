#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 12:10:56 2026

@author: bani4ki
"""
#Unit testing of minimum two components
import sys
sys.path.append("..") #To be able to import modules that reside in the parent directory
import unittest
from dictionaries import amino_acids, motifs
from functions import amino_acid_composition, calculate_properties

class TestAAComposition(unittest.TestCase):

    def test_normal_sequence(self):
        self.assertEqual(amino_acid_composition("AAAALVTVII", amino_acids),
                         {
                             'A_count': 4, 'A_fq': 0.4, 
                             'R_count': 0, 'R_fq': 0.0,
                             'N_count': 0, 'N_fq': 0.0,
                             'D_count': 0, 'D_fq': 0.0,
                             'C_count': 0, 'C_fq': 0.0,
                             'E_count': 0, 'E_fq': 0.0,
                             'Q_count': 0, 'Q_fq': 0.0,
                             'G_count': 0, 'G_fq': 0.0,
                             'H_count': 0, 'H_fq': 0.0,
                             'I_count': 2, 'I_fq': 0.2,
                             'L_count': 1, 'L_fq': 0.1,
                             'K_count': 0, 'K_fq': 0.0,
                             'M_count': 0, 'M_fq': 0.0,
                             'F_count': 0, 'F_fq': 0.0,
                             'P_count': 0, 'P_fq': 0.0,
                             'S_count': 0, 'S_fq': 0.0,
                             'T_count': 1, 'T_fq': 0.1,
                             'W_count': 0, 'W_fq': 0.0,
                             'Y_count': 0, 'Y_fq': 0.0,
                             'V_count': 2, 'V_fq': 0.2,
                             'U_count': 0, 'U_fq': 0.0 #X's count and frequency should not be included (unknown AA), 
                             }#but it should be included in the length when calculating fq of the other AAs.
                         )

    def test_only_X(self):
        self.assertEqual(amino_acid_composition("XXXXXXXXXX", amino_acids),
                         {
                             'A_count': 0, 'A_fq': 0.0, 
                             'R_count': 0, 'R_fq': 0.0,
                             'N_count': 0, 'N_fq': 0.0,
                             'D_count': 0, 'D_fq': 0.0,
                             'C_count': 0, 'C_fq': 0.0,
                             'E_count': 0, 'E_fq': 0.0,
                             'Q_count': 0, 'Q_fq': 0.0,
                             'G_count': 0, 'G_fq': 0.0,
                             'H_count': 0, 'H_fq': 0.0,
                             'I_count': 0, 'I_fq': 0.0,
                             'L_count': 0, 'L_fq': 0.0,
                             'K_count': 0, 'K_fq': 0.0,
                             'M_count': 0, 'M_fq': 0.0,
                             'F_count': 0, 'F_fq': 0.0,
                             'P_count': 0, 'P_fq': 0.0,
                             'S_count': 0, 'S_fq': 0.0,
                             'T_count': 0, 'T_fq': 0.0,
                             'W_count': 0, 'W_fq': 0.0,
                             'Y_count': 0, 'Y_fq': 0.0,
                             'V_count': 0, 'V_fq': 0.0,
                             'U_count': 0, 'U_fq': 0.0
                             }
                        )
        
    def test_new_AA(self):
        self.assertRaises(ValueError, amino_acid_composition, "AAAZ", amino_acids)
        
class TestAAProperties(unittest.TestCase):

    def test_normal_sequence(self):
        self.assertEqual(calculate_properties("AAAALVTVII", amino_acids, motifs),
                         {
                             'Hydrophobic_fq': 0.9, 'Aromatic_fq': 0.0,
                             'Charged_fq': 0.0, 'DRY_count': 0, 'DRY_fq': 0.0,
                             'NP..Y_count': 0, 'NP..Y_fq': 0.0, 'CW.P_count': 0, 
                             'CW.P_fq': 0.0, 'G.G..G_count': 0, 'G.G..G_fq': 0.0,
                             'VAIK_count': 0, 'VAIK_fq': 0.0, 'HRD_count': 0, 
                             'HRD_fq': 0.0, 'DFG_count': 0, 'DFG_fq': 0.0,
                             'APE_count': 0, 'APE_fq': 0.0, 'A.K_count': 0,
                             'A.K_fq': 0.0, 'L......L......L......L_count': 0,
                             'L......L......L......L_fq': 0.0,
                             'KIK.LK.KK_count': 0, 'KIK.LK.KK_fq': 0.0
                             }
                         )

    def test_only_X(self):
        self.assertEqual(calculate_properties("AXXXXXXXXX", amino_acids, motifs),
                         {
                             'Hydrophobic_fq': 0.1, 'Aromatic_fq': 0.0, #Handling X while including them in length of the sequence to calculate frequencies
                             'Charged_fq': 0.0, 'DRY_count': 0, 'DRY_fq': 0.0,
                             'NP..Y_count': 0, 'NP..Y_fq': 0.0, 'CW.P_count': 0,
                             'CW.P_fq': 0.0, 'G.G..G_count': 0, 'G.G..G_fq': 0.0,
                             'VAIK_count': 0, 'VAIK_fq': 0.0, 'HRD_count': 0, 
                             'HRD_fq': 0.0, 'DFG_count': 0, 'DFG_fq': 0.0, 
                             'APE_count': 0, 'APE_fq': 0.0, 'A.K_count': 0, 
                             'A.K_fq': 0.0, 'L......L......L......L_count': 0,
                             'L......L......L......L_fq': 0.0,
                             'KIK.LK.KK_count': 0, 'KIK.LK.KK_fq': 0.0
                             }
                        )

    def test_new_AA(self):
        self.assertEqual(calculate_properties("AAAZ", amino_acids, motifs),
                         {
                             'Hydrophobic_fq': 0.75, 'Aromatic_fq': 0.0, #Handling unknown AA while including them in length of the sequence to calculate frequencies
                             'Charged_fq': 0.0, 'DRY_count': 0, 'DRY_fq': 0.0,
                             'NP..Y_count': 0, 'NP..Y_fq': 0.0, 'CW.P_count': 0,
                             'CW.P_fq': 0.0, 'G.G..G_count': 0, 'G.G..G_fq': 0.0,
                             'VAIK_count': 0, 'VAIK_fq': 0.0, 'HRD_count': 0, 
                             'HRD_fq': 0.0, 'DFG_count': 0, 'DFG_fq': 0.0, 
                             'APE_count': 0, 'APE_fq': 0.0, 'A.K_count': 0, 
                             'A.K_fq': 0.0, 'L......L......L......L_count': 0,
                             'L......L......L......L_fq': 0.0,
                             'KIK.LK.KK_count': 0, 'KIK.LK.KK_fq': 0.0
                             }
                        )


if __name__ == "__main__":
    unittest.main()