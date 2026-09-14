#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 14:31:05 2026

@author: bani4ki
"""

import pandas as pd
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from functions import protein_classification, classifier
from dictionaries import amino_acids, motifs


#df = pd.read_csv("/Users/bani4ki/Documents/programming_project_bioinformatics/datasets/big_dataset.csv")
df_A= "/Users/bani4ki/Documents/programming_project_bioinformatics/datasets/dataset_A.csv"
df_B = "/Users/bani4ki/Documents/programming_project_bioinformatics/datasets/dataset_B.csv"

#df_B=df_B.drop(columns="Functional_class")

y_pred, y_test, results_df, model, label_encoder, scaler= protein_classification(df_A, amino_acids, motifs)

print(results_df.head())
print(accuracy_score(y_test, y_pred))
print(classification_report(y_test,y_pred))
print(confusion_matrix(y_test,y_pred))

# Test new dataset on an already trained model
print("new dataset:")
y_pred, y, results_df = classifier(df_B, model, amino_acids, motifs, label_encoder, scaler)
print(results_df.head())
if y is not None:
    print(accuracy_score(y, y_pred))
    print(classification_report(y,y_pred))
    print(confusion_matrix(y,y_pred)) 
else: 
    print("\n Cannot evaluate classifier performance: Unknown actual class of proteins.")

