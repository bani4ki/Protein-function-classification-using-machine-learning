#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 15:56:25 2026

@author: bani4ki
"""
#Python dictionary containing amino acids and their properties of type A : {name: int, polarity: str, charge: str, aromatic: bool, hydrophobic: bool}

amino_acids = {
    "A": {"name": "Alanine", "three_letter": "Ala", "weight": 71.04, "polarity": "nonpolar", "charge": "neutral", "aromatic": False, "hydrophobic": True},
    "R": {"name": "Arginine", "three_letter": "Arg", "weight": 156.10, "polarity": "polar", "charge": "positive", "aromatic": False, "hydrophobic": False},
    "N": {"name": "Asparagine", "three_letter": "Asn", "weight": 114.04, "polarity": "polar", "charge": "neutral", "aromatic": False, "hydrophobic": False},
    "D": {"name": "Aspartic acid", "three_letter": "Asp", "weight": 115.03, "polarity": "charged", "charge": "negative", "aromatic": False, "hydrophobic": False},
    "C": {"name": "Cysteine", "three_letter": "Cys", "weight": 103.01, "polarity": "nonpolar", "charge": "neutral", "aromatic": False, "hydrophobic": True},
    "E": {"name": "Glutamic acid", "three_letter": "Glu", "weight": 129.04, "polarity": "charged", "charge": "negative", "aromatic": False, "hydrophobic": False},
    "Q": {"name": "Glutamine", "three_letter": "Gln", "weight": 128.06, "polarity": "polar", "charge": "neutral", "aromatic": False, "hydrophobic": False},
    "G": {"name": "Glycine", "three_letter": "Gly", "weight": 57.02, "polarity": "nonpolar", "charge": "neutral", "aromatic": False, "hydrophobic": True},
    "H": {"name": "Histidine", "three_letter": "His", "weight": 137.06, "polarity": "charged", "charge": "positive", "aromatic": True, "hydrophobic": False},
    "I": {"name": "Isoleucine", "three_letter": "Ile", "weight": 113.08, "polarity": "nonpolar", "charge": "neutral", "aromatic": False, "hydrophobic": True},
    "L": {"name": "Leucine", "three_letter": "Leu", "weight": 113.08, "polarity": "nonpolar", "charge": "neutral", "aromatic": False, "hydrophobic": True},
    "K": {"name": "Lysine", "three_letter": "Lys", "weight": 128.09, "polarity": "charged", "charge": "positive", "aromatic": False, "hydrophobic": False},
    "M": {"name": "Methionine", "three_letter": "Met", "weight": 131.04, "polarity": "nonpolar", "charge": "neutral", "aromatic": False, "hydrophobic": True},
    "F": {"name": "Phenylalanine", "three_letter": "Phe", "weight": 147.07, "polarity": "nonpolar", "charge": "neutral", "aromatic": True, "hydrophobic": True},
    "P": {"name": "Proline", "three_letter": "Pro", "weight": 97.05, "polarity": "nonpolar", "charge": "neutral", "aromatic": False, "hydrophobic": True},
    "S": {"name": "Serine", "three_letter": "Ser", "weight": 87.03, "polarity": "polar", "charge": "neutral", "aromatic": False, "hydrophobic": False},
    "T": {"name": "Threonine", "three_letter": "Thr", "weight": 101.05, "polarity": "polar", "charge": "neutral", "aromatic": False, "hydrophobic": False},
    "W": {"name": "Tryptophan", "three_letter": "Trp", "weight": 186.08, "polarity": "nonpolar", "charge": "neutral", "aromatic": True, "hydrophobic": True},
    "Y": {"name": "Tyrosine", "three_letter": "Tyr", "weight": 163.06, "polarity": "polar", "charge": "neutral", "aromatic": True, "hydrophobic": False}, # Polar OH- group makes its sidechain overall hydrophilic
    "V": {"name": "Valine", "three_letter": "Val", "weight": 99.07, "polarity": "nonpolar", "charge": "neutral", "aromatic": False, "hydrophobic": True},
    "U": {"name": "Selenocysteine", "three_letter": "Sec", "weight": 150.05, "polarity": "polar", "charge": "neutral", "aromatic": False, "hydrophobic": False},
    "X": {"name": "Unknown", "three_letter": "Xaa", "weight": None, "polarity": None, "charge": None, "aromatic": None, "hydrophobic": None}
}

# Python list containing protein motifs in regex format.
motifs = ["DRY","NP..Y","CW.P","G.G..G","VAIK","HRD","DFG","APE","A.K","L......L......L......L","KIK.LK.KK"]
