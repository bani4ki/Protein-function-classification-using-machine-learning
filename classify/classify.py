#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 10:29:37 2026

@author: bani4ki
"""

import argparse
from functions import protein_classification, classifier, create_class_results_file, performance_evaluation_display, sequence_classifier
from dictionaries import amino_acids, motifs
import os



def main():
    
    parser = argparse.ArgumentParser(description = "Protein function classifier tool")
    
    #Required arguments
    parser.add_argument("dataset_path", help = "Path to protein function csv file")
    
    #Optional arguments
    parser.add_argument("--test-size",type=float, default=0.2,
                        help="Portion of data to be used as a training dataset as a float")
    parser.add_argument("--random-state-split", type=int, default=42,
                        help="Random seed of the data train/test split as an integer")
    parser.add_argument("--bootstrap", type=bool, default=False,
                        help="Bootstrap as a boolean")
    parser.add_argument("--max-depth", type=int, default=20,
                        help="Max depth of the RF classifier as an integer")
    parser.add_argument("--random-state-rf", type=int, default=42,
                        help="Random seed of the RF classifier as an integer")
    #Optional arguments - additional data to test
    parser.add_argument("--predict-class",
                        help="Path to dataset to be classified by trained model")
    parser.add_argument("--predict-sequence",
                        help="Protein amino acid sequence to be classified by trained model as a string")
    
    #Parse agument from command line
    args = parser.parse_args()
    try:
        results_df, model, label_encoder, scaler= protein_classification(args.dataset_path, 
                                                                         amino_acids,
                                                                         motifs,
                                                                         args.test_size, 
                                                                         args.random_state_split, 
                                                                         args.bootstrap, 
                                                                         args.max_depth,
                                                                         args.random_state_rf)
    except FileNotFoundError:
        print(f"Error: The file '{args.dataset_path}' could not be found. \n Please check the file name and try again.")
        return
    except KeyError as e:
        print(f"Error: dataset is missing a required column: {e}")
        return
    
    if args.predict_class:
            
        print("\n Known dataset classification results: \n ")
    else:
        print("\n Dataset classification results: \n ")
                
    performance_evaluation_display(results_df, label_encoder)
            
    results_dir = "results"
    train_test_dir = os.path.join(results_dir,"train_test_results")
    os.makedirs(train_test_dir, exist_ok=True)
    create_class_results_file(results_df, output_directory= train_test_dir)
            
    
    
         
    if args.predict_class:
        
        try:
            results_df = classifier(args.predict_class, model, amino_acids, motifs, label_encoder, scaler)
            print("\n Unknown dataset classification results: \n ")
            
            performance_evaluation_display(results_df, label_encoder)
            
            implementation_dir = os.path.join(results_dir,"implementation_results")
            os.makedirs(implementation_dir, exist_ok=True)
                    
            create_class_results_file(results_df, output_directory=implementation_dir)
            
        except FileNotFoundError:
            print(f"Error: The file '{args.predict_class}' could not be found. \n Please check the file name and try again.")
        except KeyError as e:
            print(f"Error: dataset is missing a required column: {e}")
            return
        except Exception as e:
            print(f"Unexpected internal error while classifying dataset: {e}")  #internal error handling
    if args.predict_sequence:
        try:
            y_pred = sequence_classifier(args.predict_sequence, model, amino_acids, motifs, label_encoder, scaler)
            print(f"\n Predicted functional class of the sequence: {y_pred}")
        except ValueError as e:
            print(f"Error: {e}")
            return   
        
if __name__ == "__main__":
    main()