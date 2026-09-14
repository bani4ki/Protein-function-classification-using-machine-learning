#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 18:55:49 2026

@author: bani4ki
"""

import pandas as pd
from sklearn.model_selection import train_test_split

df = pd.read_csv("/Users/bani4ki/Documents/programming_project_bioinformatics/datasets/uniprotkb_keyword_KW_0418_OR_keyword_KW_2026_09_07.tsv", sep="\t")

keywords = ["Kinase", "G protein-coupled receptor",
            "Transcription", "Ribosomal protein"]
   
def find_functional_class(keyword_string):
    if pd.isna(keyword_string):
        return None
    for keyword in keywords:
        if keyword in keyword_string:
            return keyword
    return None

df['Functional_class'] = df['Keywords'].apply(find_functional_class)

clean_df = df.dropna(subset=['Functional_class', 'Sequence'])

clean_df[['Entry','Sequence', 'Functional_class']].to_csv("big_dataset.csv", index=False)


df_clean = pd.read_csv("big_dataset.csv")

dataset_A, dataset_B = train_test_split(
    df_clean, 
    test_size=0.5, 
    random_state=42, 
    stratify=df_clean["Functional_class"]
)

_, dataset_dummy = train_test_split(
    df_clean, 
    test_size=0.01, 
    random_state=42, 
    stratify=df_clean["Functional_class"]
)

dataset_A.to_csv("dataset_A.csv", index=False)
dataset_B.to_csv("dataset_B.csv", index=False)
dataset_dummy.to_csv("dataset_dummy.csv", index=False)

if __name__ == "__main__":
    
    print("Cleaned data functional classes separation:")
    print(clean_df['Functional_class'].value_counts(), "\n")
    
    print(f"Rows and columns of dataset A: {dataset_A.shape}")
    print(f"Rows and columns of dataset B: {dataset_B.shape}")
    print(f"Rows and columns of dummy dataset: {dataset_dummy.shape}\n")
    print("Dummy dataset functional classes separation:")
    print(dataset_dummy['Functional_class'].value_counts())
print(df.head())