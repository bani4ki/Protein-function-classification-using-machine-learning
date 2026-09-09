#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 14:31:27 2026

@author: bani4ki
"""
from dictionaries import amino_acids
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.ensemble import RandomForestClassifier


def amino_acid_composition(sequence):
    metrics = {}
    sequence = sequence.upper()
    
    for aa in amino_acids.keys():
        count = sequence.count(aa)
        metrics[f"{aa}_count"] = count
        frequency = count/len(sequence)
        metrics[f"{aa}_fq"] = round(frequency, 2)    

    return metrics

def calculate_properties(sequence):
    metrics = {}
    sequence = sequence.upper()

    hydrophobic = sum(
        1 for aa in sequence if aa in amino_acids and amino_acids[aa].get("hydrophobic") == True
    )
    metrics["Hydrophobic_fq"] = round(hydrophobic/len(sequence),2)
    
    aromatic = sum(
        1 for aa in sequence if aa in amino_acids and amino_acids[aa].get("aromatic") == True
    )
    metrics["Aromatic_fq"] = round(aromatic/len(sequence),2)
    
    charged = sum(
        1 for aa in sequence if aa in amino_acids and amino_acids[aa].get("charge") != "neutral"
    )
    metrics["Charged_fq"] = round(charged/len(sequence),2)
    
    return metrics

def dataset_preparation_for_rf(df):
    df["Sequence_length"] = df["Sequence"].str.len()
    amino_acid_props_df = df["Sequence"].apply(calculate_properties).apply(pd.Series)
    amino_acid_comp_df = df["Sequence"].apply(amino_acid_composition).apply(pd.Series)
    df = pd.concat([df, amino_acid_comp_df, amino_acid_props_df], axis = 1)
    df_rf = df.loc[ : , (df.columns != "Sequence")]
    return df_rf

def protein_classification(dataset, test_size=0.2, random_state_split=42, bootstrap=False, max_depth=20,random_state_rf = 42):
    
    df_rf = dataset_preparation_for_rf(dataset)
    X,y = df_rf.loc[ : , (df_rf.columns != "Functional_class")], df_rf["Functional_class"]
    X_train, X_test, y_train, y_test = train_test_split(X,y,test_size = test_size, random_state = random_state_split)
   
    label_encoder = LabelEncoder()
    y_train = label_encoder.fit_transform(y_train)
    y_test = label_encoder.transform(y_test)
    
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train) 
    X_test = scaler.transform(X_test)
    
    rf_model = RandomForestClassifier(bootstrap=bootstrap, max_depth=max_depth,random_state = random_state_rf)
    rf_model.fit(X_train,y_train)
    y_pred = rf_model.predict(X_test)
    
    mapping_df = df_rf[['Functional_class']].copy() 
    mapping_df['Functional_class_encoded'] = label_encoder.transform(df_rf['Functional_class'].values)
    mapping_df = mapping_df.drop_duplicates('Functional_class')
    mapping_df = mapping_df.sort_values('Functional_class_encoded',ignore_index = True)
    
    return y_pred, y_test, mapping_df



