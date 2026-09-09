#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 14:31:05 2026

@author: bani4ki
"""

import pandas as pd
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from functions import protein_classification


df = pd.read_csv("/Users/bani4ki/Documents/programming_project_bioinformatics/datasets/big_dataset.csv")
#df = pd.read_csv("/Users/bani4ki/Documents/programming_project_bioinformatics/datasets/dataset_A.csv")
#df = pd.read_csv("/Users/bani4ki/Documents/programming_project_bioinformatics/datasets/dataset_B.csv")


y_pred, y_test, mapping_df = protein_classification(df)

print(mapping_df)
print(accuracy_score(y_test, y_pred))
print(classification_report(y_test,y_pred))
print(confusion_matrix(y_test,y_pred))

