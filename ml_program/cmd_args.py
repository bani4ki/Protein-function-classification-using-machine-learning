#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 10:29:37 2026

@author: bani4ki
"""

import argparse
from functions import protein_classification, classifier, create_class_results_file, performance_evaluation_display, sequence_classifier
from dictionaries import amino_acids, motifs


def main():
    
    parser = argparse.ArgumentParser(description = "Protein function classifier tool")
    
    #Required arguments
    parser.add_argument("dataset_path", help = "Path to protein function csv file")
    
    #Optional arguments
    parser.add_argument("--test-size", default=0.2,
                        help="Portion of data to be used as a training dataset as a float")
    parser.add_argument("--random-state-split", default=42,
                        help="Random seed of the data train/test split as an integer")
    parser.add_argument("--bootstrap", default=False,
                        help="Bootstrap as a boolean")
    parser.add_argument("--max-depth", default=20,
                        help="Max depth of the RF classifier as an integer")
    parser.add_argument("--random-state-rf", default=42,
                        help="Random seed of the RF classifier")
    #Optional arguments - additional data to test
    parser.add_argument("--dataset-to-classify",
                        help="Path to dataset to be classified by trained model")
    parser.add_argument("--protein-to-classify",
                        help="Protein amino acid sequence to be classified by trained model as a string")
    
    #Parse agument from command line
    args = parser.parse_args()
    
    results_df, model, label_encoder, scaler= protein_classification(args.dataset_path, 
                                                                     amino_acids,
                                                                     motifs,
                                                                     args.test_size, 
                                                                     args.random_state_split, 
                                                                     args.bootstrap, 
                                                                     args.max_depth,
                                                                     args.random_state_rf)
    performance_evaluation_display(results_df, label_encoder)
    
    create_class_results_file(results_df)
    
    if args.dataset_to_classify:
        results_df = classifier(args.dataset_to_classify, model, amino_acids, motifs, label_encoder, scaler)
    
        performance_evaluation_display(results_df, label_encoder)
    
        create_class_results_file(results_df)
    
    if args.protein_to_classify:
        y_pred = sequence_classifier(args.protein_to_classify, model, amino_acids, motifs, label_encoder, scaler)
        print(y_pred)
    
if __name__ == "__main__":
    main()