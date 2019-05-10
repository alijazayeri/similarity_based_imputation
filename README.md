# Authors: Ali Jazayeri, Ou Stella Liang, and Christopher C. Yang
# Date: 05/10/2019

# This code is developed for the IEEE ICHI - 2019 Data Analytics Challenge on Missing data Imputation (DACMI).
# Details are provided at http://www.ieee-ichi.org/challenge.html

# The code for the Similarity-based Imputation of Missing Data in Electronic Health Records is written in R. 
# To run the code, the folders of training and testing data sets should be added to the same directory that the R file is placed, with the name of:
#  1. "training_data", a folder containing the training data files (reference files)
#  2. "testing_data", a folder containing the test data files
# Also, the coefficients.csv file should be placed in the R file directory. This file contains the parameters of generalized bell-shaped functions used in the code.
# Then, the code creates a new folder, "imputed_results", for .
# The code creates one single data frame composed of all the reference files.
# This data frame is used as a reference in the code for imputing missing values in test files.
# After all the missing data in each test file are imputed, the test file is exported to the "imputed_results" with the same name it has had in the "testing_data" folder
























