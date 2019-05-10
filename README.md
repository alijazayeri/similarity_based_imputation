# Authors: Ali Jazayeri, Ou Stella Liang, and Christopher C. Yang
# Date: 05/10/2019

# This code is developed for the IEEE ICHI - 2019 Data Analytics Challenge on Missing data Imputation (DACMI).
# Details are provided at http://www.ieee-ichi.org/challenge.html

# The code for the Similarity-based Imputation of Missing Data in Electronic Health Records is written in R. 
# To run the code, you need two folders in the same directory that the code file is placed:
#  1. "training_data", a folder containing the training data (reference files)
#  2. "testing_data", a folder containing the test data
# Then, one new folder, "imputed_results", is automatically created.
# Also the coefficients.csv file should be placed in the code file directory. This file contains the generalized bell-shaped functions used in the code.
# The code creates one single data frame composed of all the reference files.
# This data frame is used as a reference in the code for imputing missing values in test files.
# After all the missing data in each test file are imputed, the test file is exported to the "imputed_results" with the same name it has had in the "testing_data" folder
























