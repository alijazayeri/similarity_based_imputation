# importing some libraries
library(dplyr)
library(ggplot2)
library(pdist)  # --> this one will be used for computing the distance between two vectors
library(zoo)    # --> this one will be used for linear interpolation of missing values
library(stringi)


# Getting and setting the working directory, and the paths to the training and testing data
current_directory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_directory)
training_data_directory <- paste(current_directory,'/training_data',sep = '')
testing_data_directory <- paste(current_directory,'/testing_data',sep = '')

dir.create(file.path(current_directory,'imputed_results')) # --> Creates an empty folder for imputed test results

#_______________________________________________________________________________
# Creating (and exporting) one single data frame composed of all patients' data
# For each patient an "id" will be added as a new column
# This block needs to be executed just once, the generated data frame will be exported
# to the current directory 
# It is commented out now!
# traing_files_list <- list.files(path = training_data_directory, pattern = '.csv', full.names = TRUE)
# 
# read_training_data <- function(x){
#   temp_tr_df <- read.table(x, header = TRUE,sep = ',')
#   temp_id <- as.numeric(stri_match_last_regex(x,'\\d+'))
#   temp_tr_df$id <- temp_id
#   return(temp_tr_df)
# }
# 
# train_df <- do.call('rbind', lapply(traing_files_list, function(x) {
#   read_training_data(x)
# }))
# train_df <- train_df[order(train_df$id),]
# 
# write.csv(train_df,'train_dataframe.txt',row.names = FALSE)
# End of creating the "training" data frame
#_______________________________________________________________________________

# Importing the single data frame composed of the data of all the patients in the training dataset
# This data frame will be used as a reference for patients in the testing dataset
train_df <- read.csv('train_dataframe.txt')

# Columns need to be saved
orig_cols <- names(train_df)[-15]

# List of analytes (lab tests)
analyte_list <- setdiff(names(train_df),c("CHARTTIME","id"))

# Add two columns for computing the previous and next rows' time difference with the current row
train_df <- train_df %>% group_by(id) %>% 
  mutate(dt_mOne = 
           ifelse(row_number() == 1,CHARTTIME[2]  , ifelse(row_number() == n(),CHARTTIME[n()-1]-CHARTTIME[n()-2], 
                                                           CHARTTIME - lag(CHARTTIME))),
         dt_pOne =  ifelse(row_number() == 1,CHARTTIME[3]-CHARTTIME[2],ifelse(row_number() == n(),CHARTTIME[n()]-CHARTTIME[n()-1],lead(CHARTTIME) - CHARTTIME))) %>% 
  ungroup

# Importing the coefficients for the similarity calculation
coef_df <- read.csv('coefficients.csv')

# Here, we just include the rows for comparison which have more than 9 (out of 13) known values for comparison
# the "9" is arbitrary, and we didn't evaluate other values
# Basically, it is assumed that rows having less than 9 known values are not reliable enough to be used for
# estimating missed values
train_df$missed_df_non_na <- apply(train_df[,analyte_list], 1, function(x) length(x[!is.na(x)]))
train_df$missed_df_non_na <- ifelse(train_df$missed_df_non_na > 9, 1,NA)

# This function computes the weights based on the distance
# Basically, this is a generalized bell-shaped function with three parameters a,b,c
# which gets a vector and computes a weight for each element in the vector based on its value
weight_func <- function(ddt,a,b,c){
  return(sapply(ddt, function(x) ifelse(is.na(x),100000,1/(1+(abs((x-c)/a)^(2*b))))))
}

# List of test dataset files
testing_files_list <- list.files(path = testing_data_directory, pattern = '.csv', full.names = TRUE)

for (test_file in testing_files_list[3000:3999]){
  temp_test_df <- read.table(test_file, header = TRUE,sep = ',')
  temp_id <- as.numeric(stri_match_last_regex(test_file,'\\d+'))
  temp_test_df$id <- temp_id

  # Two functions created to compute a range for each analyte of the patient
  # These functions help to evaluate the acceptance of each estimated missing value
  iqr_fun_max <- function(x) {max(x, na.rm=TRUE) + 0.25 * IQR(x, na.rm=TRUE)}
  iqr_fun_min <- function(x) {min(x, na.rm=TRUE) - 0.25 * IQR(x, na.rm=TRUE)}
  
  # Creating two data frames providing the acceptance range for estimated values 
  iqr_df_max <- as.data.frame(t(apply(temp_test_df[,analyte_list], 2, iqr_fun_max)))
  iqr_df_min <- as.data.frame(t(apply(temp_test_df[,analyte_list], 2, iqr_fun_min)))
  
  # Adding two columns to the patient data frame which shows the delta-t 
  # or time difference for each row relative to the previous and next measurement
  temp_test_df$dt_mOne <- temp_test_df$CHARTTIME - lag(temp_test_df$CHARTTIME)
  temp_test_df$dt_mOne[1] <- temp_test_df$CHARTTIME[2]
  temp_test_df$dt_pOne[1] <- temp_test_df$CHARTTIME[3] - temp_test_df$CHARTTIME[2]
  temp_test_df$dt_pOne <- lead(temp_test_df$CHARTTIME) - temp_test_df$CHARTTIME
  temp_test_df$dt_mOne[nrow(temp_test_df)] <- temp_test_df$CHARTTIME[nrow(temp_test_df)-1] - temp_test_df$CHARTTIME[nrow(temp_test_df)-2]
  temp_test_df$dt_pOne[nrow(temp_test_df)] <- temp_test_df$CHARTTIME[nrow(temp_test_df)] - temp_test_df$CHARTTIME[nrow(temp_test_df)-1]
  # We impute missing values row-wise. 
  temp_test_df$na_cnt <- apply(temp_test_df[,analyte_list], 1, function(x) length(which(is.na(x))))
  temp_test_df$row_n <- seq(1:nrow(temp_test_df))
  # The following data frame contains the row number of the patient data frame with the missing analytes
  # in the row. Just rows with missing values are included. Rows with multiple missing values are
  # repeated with the same row number if row_n column and different analyte names in analyte column
  temp_test_na_df <- bind_rows(apply(temp_test_df[which(temp_test_df$na_cnt>0),],1, function(x) 
    data.frame(row_n = x[length(x)], analyte = unlist(intersect(names(x[which(is.na(x))]),analyte_list)))))
  
  # The function impute_values called with the reference data frame, patient data frame, and missing values data frame
  temp_imputed_df <- impute_values(train_df, temp_test_df, temp_test_na_df)
  
  # Export the patient data frame with imputed missing values to the "imputed_results" folder
  write.csv(temp_imputed_df[,orig_cols],paste(current_directory,'/imputed_results/',temp_id,'.csv',sep = ''),row.names = FALSE)
}

#____________________________________________________________________________________________________
# This function is the main contribution of this work!
# It estimates the missing data based on other patients' data, if there is not enough-similar patient
# in the reference data frame (training data), it simply interpolates the missing data linearly!
impute_values <- function(ref_df, pat_df, miss_df){
  
  # Iteration over rows with missing values
  rows_w_missing <- unique(miss_df$row_n)
  for (r in rows_w_missing){
    # Computing the distance between the row of the patient with all the rows in the reference data frame
    # Note 1: if the row in the patient data frame is the first row, we consider the similarity of the first three row
    #       if the row in the patient data frame is the last row, we consider the similarity of the last three row
    if (r == 1){
      r_mOne <- r + 1
      r_pOne <- r + 2
    }
    else if (r == nrow(pat_df)){
      r_mOne <- r - 2
      r_pOne <- r - 1      
    }
    else{
      r_mOne <- r - 1
      r_pOne <- r + 1
    }
    r_min <- max(r-2,0)
    r_max <- min(r+2,nrow(pat_df))
    
    
    # The distances between the row with missing value and all the rows in the reference data frame
    # is computed. Also, the similarity of the previous and next rows are also computed,
    # (unless, it is the first or last row, in this case, refer to the Note 1 above)
    dist <- pdist(pat_df[r,analyte_list],train_df[,analyte_list])@dist
    dist_mOne <- pdist(pat_df[r_mOne,analyte_list],train_df[,analyte_list])@dist
    dist_pOne <- pdist(pat_df[r_pOne,analyte_list],train_df[,analyte_list])@dist
    ddt_mOne <- pat_df$dt_mOne[r] - train_df$dt_mOne
    ddt_pOne <- pat_df$dt_pOne[r] - train_df$dt_pOne
    
    # Iterating over analytes
    for (a in miss_df$analyte[which(miss_df$row_n==r)]){
      # The parameters of generalized bell-shaped function are retrieved based on the 
      # analyte with missing value in the patient data frame
      wd <- coef_df$a_row[which(coef_df$analyte==a)]
      sp <- coef_df$b_row[which(coef_df$analyte==a)]
      wd_dist <- coef_df$a_dist[which(coef_df$analyte==a)]
      sp_dist <- coef_df$b_dist[which(coef_df$analyte==a)]
      c <- 0
      
      # The weight of different rows are computed based on the differences between
      # inter-measurement times
      wt_mOne <- weight_func(ddt_mOne,wd,sp,c)
      wt_pOne <- weight_func(ddt_pOne,wd,sp,c)
      
      # The total distance is computed based on the status of the patient in row r
      # and previous and next rows (unless, it is the first or last row, in this case, refer to the Note 1 above)
      dist_tot <- dist + dist_mOne * wt_mOne + dist_pOne * wt_pOne
      
      # The minimum distance is considered for the c parameter of the bell-shaped function
      # for computing the estimated missing value
      minimum_dist <- min(dist_tot)
      dist_wt <- weight_func(dist_tot,wd_dist,sp_dist,minimum_dist)
      
      # Here, we just include rows with known values for the desired analyte in the reference data frame
      # and also just the rows with more than 9 known analyte values ("9" explained above)
      na_idx <- which(!is.na(train_df[,a]) & !is.na(train_df$missed_df_non_na))
      
      # The estimated value for missing data point is calculated as a weighting average of 
      # the known values in the patients of the reference data frame based on their similarity
      # to the current patient
      weighted_avg <- sum(train_df[na_idx,a]*dist_wt[na_idx]/sum(dist_wt[na_idx]), na.rm = TRUE)
      
      
      # If the estimated value is significantly out of the range of measured values for the
      # desired analyte of the patient, it is assumed that enough-similar patient is not found,
      # then, we just use linear interpolation to estimate the missing value, simple, right? :)
      iqr_max <- max(pat_df[r_min:r_max,a], na.rm=TRUE) + 0.25 * IQR(pat_df[r_min:r_max,a], na.rm=TRUE)
      iqr_min <- min(pat_df[r_min:r_max,a], na.rm=TRUE) - 0.25 * IQR(pat_df[r_min:r_max,a], na.rm=TRUE)
      
      # in cases where there is not enough neighbor (all the neighbors are NA), we just
      # consider iqr_max and iqr_min to be weighted_avg to prevent from any error
      iqr_max <- ifelse(!is.na(iqr_max),iqr_max, weighted_avg)
      iqr_min <- ifelse(!is.na(iqr_min),iqr_min, weighted_avg)
      
      if (weighted_avg > iqr_max |
          weighted_avg < iqr_min){
        
        
        # For linear interpolation we use na.approx from "zoo" library
        weighted_avg <- na.approx(pat_df[,a],
                                  pat_df[,c('CHARTTIME')],na.rm = FALSE,rule=2)[r]

      }
      
      # The value of missing data points is placed in the patient data frame
      pat_df[r,a] <- weighted_avg
    }
  }
  return(pat_df)
}
#____________________________________________________________________________________________________