# Protein-function-classification-using-machine-learning
A command line tool to classify protein into functional classes using machine learning (Random Forest Algorithm).  
## Author
Yana Atanasova     
## Description
This program will take a dataset of protein sequences and their functions as input, split it in a proportion provided by the user, train, test, and evaluate a Random Forest machine learning model with adjustable parameters on the dataset , and classify and evaluate additional datasets and sequences if provided with such.    

---   

## Installation
Clone the repository and navigate to the directory:   
```bash
git clone git@github.com:bani4ki/Protein-function-classification-using-machine-learning.git
cd classify
```   
### Dependencies
- Python version: Python 3.9 or later (developed and tested on 3.13)   
  Check Python version:   
  ```
    python --version
    ```   
    
    And if needed, update via the official latest installer from Python Downloads.   
  
- Libraries: pandas, numpy, scikit-learn      
  - Installation:   
    ```
    pip install pandas numpy scikit-learn
    ```   
## Usage
Run from the command line, providing the required parameter:      
```bash
python classify.py path/to/dataset    
```
- Required parameter: path to a training dataset in csv format, e.g. a Uniprot dataset including a Functional_class column. These datasets can be obtained by following instructions in the [supplementary materials](./classify/supplementary_materials) folder.    

Or, run using optional parameters:     
```bash
python classify.py path/to/training/dataset --test-size 0.2 --predict-class path/to/dataset
```    
You can see all optional parameters by running `python classify.py --help`.      
- Optional parameters: Values for RF classifier parameters (defaults obtained after hyperparameter tuning)
  - `--test-size` : Portion of data to be used as a training dataset as a float. Default 0.2.
  - `--random-state-split`: Random seed of the data train/test split as an integer. Default 42.
  - `--bootstrap`: Bootstrap as a boolean. Default False.
  - `--max-depth`: Max depth of the RF classifier as an integer. Default 20
                        help="Max depth of the RF classifier as an integer")
  - `--random-state-rf`: Random seed of the RF classifier as an integer. Default 42.
- Optional parameters: Additional data to be classified or tested 
  - `--predict-class`: Path to dataset to be classified by trained model. These datasets can be obtained by following instructions in the [supplementary materials](./classify/supplementary_materials) folder.
  - `--predict-sequence`: Protein amino acid sequence to be classified by trained model as a string.
### Output 
The program will save the results to disk. It will create a directory with results from the training and testing data, and from the additional data as csv files in separate directories. It will also print out parameters of the model, performance evaluation, and a confirmation message to the console. 

---   

## Documentation
The program documentation is provided within the program files as docstrings. The hyperparameter tuning results and scripts are provided in the [supplementary materials](./classify/supplementary_materials).
## Testing
A testing documentation and results are provided in the [supplementary materials](./classify/supplementary_materials).


