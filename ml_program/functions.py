#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 14:31:27 2026

@author: bani4ki
"""
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.ensemble import RandomForestClassifier
import re
import os


def amino_acid_composition(sequence, amino_acids):
    """
    Calculates the amino acid composition and frequency 
    of the different amino acids in an amino acid sequence.

    Parameters
    ----------
    sequence : str
        A sequence of amino acids using the standard single-letter code.
    amino_acids : dictionary
        Python dictionary containing amino acids and their properties.
        
    Returns
    -------
    dictionary
        Python dictionary containing amino acid composition and frequency
        in the sequence.
    """
    metrics = {}
    sequence = sequence.upper()
    
    for aa in amino_acids.keys():
        count = sequence.count(aa)
        metrics[f"{aa}_count"] = count
        frequency = count/len(sequence)
        metrics[f"{aa}_fq"] = round(frequency, 2)    

    return metrics

def calculate_properties(sequence, amino_acids, motifs):
    """
    Calculates the properties of amino acids in a sequence by count and by frequency.

    Parameters
    ----------
    sequence : str
        A sequence of amino acids using the standard single-letter code.
    amino_acids : dictionary
        Python dictionary containing amino acids and their properties.
    motifs : list[str]
        Python list containing protein motifs in regex format.
        
    Returns
    -------
    dictionary
        Python dictionary containing amino acid properties and their frequency,
        and motif count and frequency in the sequence.
    """
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
    
    for motif in motifs:
        count = len(re.findall(motif,sequence))
        metrics[f"{motif}_count"] = count
        frequency = count/len(sequence)
        metrics[f"{motif}_fq"] = round(frequency, 2)
    
    return metrics

def dataset_preparation_for_rf(dataset, amino_acids, motifs):
    """
    Prepares a dataset or a single sequence to be classified by the RF model 
    by formatting it in a suitable way. Calculates the amino acid composition 
    and properties of the sequence(s).

    Parameters
    ----------
    dataset : OS path as str or Pandas Dataframe
        Set of identifiers, amino acid sequences and their functional class 
        in separate columns. 
    amino_acids : dictionary
        Python dictionary containing amino acids and their properties.
    motifs : list[str]
        Python list containing protein motifs in regex format.

    Returns
    -------
    Pandas DataFrame
        A dataframe with separate columns for identifiers, amino acid sequences and their 
        functional class, and their amino acid properties and composition.
    """
    if isinstance(dataset, pd.DataFrame):
        df = dataset
    else: 
        df = pd.read_csv(dataset)
    df["Sequence_length"] = df["Sequence"].str.len()
    amino_acid_props_df = df["Sequence"].apply(calculate_properties, args=(amino_acids,motifs,)).apply(pd.Series)
    amino_acid_comp_df = df["Sequence"].apply(amino_acid_composition, args=(amino_acids,)).apply(pd.Series)
    df = pd.concat([df, amino_acid_comp_df, amino_acid_props_df], axis = 1)
    df_rf = df.loc[ : , (df.columns != "Sequence")]
    return df_rf

def protein_classification(dataset, amino_acids, motifs, 
                           test_size=0.2, random_state_split=42, 
                           bootstrap=False, max_depth=20,random_state_rf = 42):
    """
    Trains and tests a Random Forest classifier to predict proteins' 
    functional class. Optional parameters are set as the calculated 
    optimal parameters after hyperparameter tuning of the model.

    Parameters
    ----------
    dataset : Pandas DataFrame
        A dataframe with separate columns for identifiers, amino acid sequences
        and functional class, and their amino acid properties and composition.
    amino_acids : dictionary
        Python dictionary containing amino acids and their properties.
    motifs : list[str]
        Python list containing protein motifs in regex format.
    test_size : float, optional
        Portion of data to be used as a training dataset. The default is 0.2.
    random_state_split : int, optional
        Random seed of the data train/test split. The default is 42.
    bootstrap : bool, optional
        Is the model using bootstrapping. The default is False.
    max_depth : int, optional
        Max depth of the RF classifier. The default is 20.
    random_state_rf : int, optional
        Random seed of the RF classifier. The default is 42.

    Returns
    -------
    results_df : Pandas DataFrame
        A dataframe with separate columns for identifiers of the test dataset 
        and their predicted functional class via the RF model, and their actual 
        functional class.
    rf_model : RandomForestClassifier
        The trained random forest model instance.
    label_encoder : LabelEncoder
        The fitted target label encoder.
    scaler : StandardScaler
        The fitted feature scaler.
    """
    
    df_rf = dataset_preparation_for_rf(dataset, amino_acids, motifs)
    X,y = df_rf.loc[ : , (df_rf.columns != "Functional_class") 
                    & (df_rf.columns != "Entry")], df_rf["Functional_class"]
    entries = df_rf["Entry"]
    X_train, X_test, y_train, y_test, entry_train, entry_test = train_test_split(
        X,y,entries,test_size = test_size, random_state = random_state_split)
    
    
    label_encoder = LabelEncoder()
    y_train = label_encoder.fit_transform(y_train)
    y_test = label_encoder.transform(y_test)

    
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train) 
    X_test = scaler.transform(X_test)
    
    rf_model = RandomForestClassifier(bootstrap=bootstrap, max_depth=max_depth,
                                      random_state = random_state_rf)
    rf_model.fit(X_train,y_train)
    y_pred = rf_model.predict(X_test)
    
    results_df = pd.DataFrame({"Entry": entry_test,
                               "Predicted_class": label_encoder.inverse_transform(y_pred),
                               "Actual_class": label_encoder.inverse_transform(y_test)})
    return results_df, rf_model,label_encoder, scaler


def classifier(dataset, model, amino_acids, motifs, label_encoder, scaler):
    """
    Fits a trained Random Forest classifier to a dataset to predict proteins' 
    functional class.

    Parameters
    ----------
    dataset : Pandas DataFrame
        A dataframe with separate columns for identifiers, amino acid sequences
        and functional class, and their amino acid properties and composition.
    model : RandomForestClassifier
        A trained random forest model.
    amino_acids : dictionary
        Python dictionary containing amino acids and their properties.
    motifs : list[str]
        Python list containing protein motifs in regex format.
    label_encoder : LabelEncoder
        A target label encoder.
    scaler : StandardScaler
        A feature scaler.

    Returns
    -------
    Pandas DataFrame
        A dataframe with separate columns for identifiers of the dataset 
        and their sequences' predicted functional class via the RF model. If 
        their actual functional class exists, it is included.
    """
    
    if dataset is not None:
        df_rf = dataset_preparation_for_rf(dataset, amino_acids, motifs)
        X = df_rf.loc[ : , (df_rf.columns != "Functional_class")
                      & (df_rf.columns != "Entry")]
        if "Functional_class" in df_rf.columns:
            y = df_rf["Functional_class"]
            y = label_encoder.transform(y)
        else: 
            y = None
        
        entries = df_rf["Entry"]
    
        X = scaler.transform(X)
        y_pred = model.predict(X)
    
    
        results_df = pd.DataFrame({"Entry": entries,
                                   "Predicted_class": label_encoder.inverse_transform(y_pred)})
        if y is not None:
            results_df["Actual_class"] = label_encoder.inverse_transform(y)
        else:
            results_df["Actual_class"] = y

    
        return results_df
    
def sequence_classifier(sequence, model, amino_acids, motifs, label_encoder, scaler):
    """
    Fits a trained Random Forest classifier to a dataset to predict proteins' 
    functional class.

    Parameters
    ----------
    sequence : str
        A sequence of amino acids using the standard single-letter code.
    model : RandomForestClassifier
        A trained random forest model.
    amino_acids : dictionary
        Python dictionary containing amino acids and their properties.
    motifs : list[str]
        Python list containing protein motifs in regex format.
    label_encoder : LabelEncoder
        A target label encoder.
    scaler : StandardScaler
        A feature scaler.

    Returns
    -------
    str
        The sequence's predicted functional class via the RF model.
    """
    df=pd.DataFrame({"Sequence": [sequence.upper()]})
    df_rf=dataset_preparation_for_rf(df, amino_acids, motifs)
    X = df_rf.loc[ : , (df_rf.columns != "Functional_class")]
    X = scaler.transform(X)
    y_pred = model.predict(X)
    return label_encoder.inverse_transform(y_pred)

def create_class_results_file(results_df):
    """
    Writes results of model classification to a csv file with unique name.

    Parameters
    ----------
    results_df : Pandas DataFrame
        A dataframe with separate columns for identifiers of a dataset,
        their predicted functional class, and their actual functional class.

    Returns
    -------
    int
        Suffix of the newly created file name.
    """
    filename = "protein_classification_results.csv"
    basename = "protein_classification_results"
    extension=".csv"
    suffix = 0
    
    while os.path.exists(filename):
        suffix += 1
        filename = f"{basename}_{suffix}{extension}"
        
    results_df.to_csv(filename, index=False)

    print(f"Results successfully written to disk in file {filename}")
    return suffix

def performance_evaluation_display(results_df, label_encoder):
    """
    Evaluates performance of a classifier by calculating its accuracy and 
    creating a classification report and confusion matrix. Prints results.

    Parameters
    ----------
    results_df : Pandas DataFrame
        A dataframe with separate columns for identifiers of a dataset,
        their predicted functional class, and their actual functional class.
    label_encoder : LabelEncoder
        A target label encoder.

    Returns
    -------
    None
    """
    y_pred=label_encoder.transform(results_df["Predicted_class"])
    y = label_encoder.transform(results_df["Actual_class"])
    print(f"Results head: \n {results_df.head()} \n ")
    if y is not None:
        print(f"Accuracy score of classifier: {round(accuracy_score(y, y_pred),2)} \n ")
        print(f"Classification report: \n {classification_report(y,y_pred)} \n ")
        print(f"Confusion matrix: \n {confusion_matrix(y,y_pred)} \n ") 
    else: 
        print("\n Cannot evaluate classifier performance: Unknown actual class of proteins.")
    
    
    