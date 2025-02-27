### This is the R-script to combine ABM output files into a single big .csv file for data analysis ###
## Note: It is not recommended to upload the .csv file > 100MB to GitHub ##

## Prepare the environment

# set working directory (for YJ's laptop: replace YingJie with lanexran)
setwd("/Users/YingJie/Documents/GitHub//ABM_code/analysis")

# packages & customized functions for data processing
source("/Users/YingJie/Documents/GitHub/ABM_code/analysis/process_util.R") 
source("/Users/YingJie/Documents/GitHub/ABM_code/analysis/figure_themes.R")

## Load the design

# Read your design CSV file into R
file_path <- "/Users/YingJie/Documents/GitHub/ABM_code/input/3x4_test_IP_ctl_K_LB6-2_traits-fixed.csv"
design <- read.csv(file_path, header = TRUE)

## Create a list to store the SQLite output
results_list <- list()

# list the treatment levels
treatment_I <- unique(design$I)
treatment_P <- unique(design$P)

# assign trt I & P 
for (j in first_key:last_key) {
  row_index = which(design$key == j)
  exp = design$exp[row_index]
  
  trtI <- ifelse(design$I[row_index] == treatment_I[1], "full",
                 ifelse(design$I[row_index] == treatment_I[2], "nested", "modular"))
  
  trtP <- ifelse(design$P[row_index] == treatment_P[1], "full",
                 ifelse(design$P[row_index] == treatment_P[2], "modular", "hub"))
  
  rep = design$rep[row_index]
  
  # Generate the .sqlite file path based on the key
  path <- sprintf("/Users/YingJie/Documents/GitHub/ABM_code/output/output%d.sqlite", j)
  
  
  # Connect to the SQLite database
  db <- dbConnect(dbDriver("SQLite"), path)
  
  # Execute the query and store the result in the list; also add columns of experiment, treatment, and replicate
  result <- dbGetQuery(db, "SELECT * FROM bsubabundance")
  n = length(result$t)
  result <- c(list(experiment = rep(exp, n),
                   group_I = rep(trtI, n),
                   group_P = rep(trtP, n),
                   replicate = rep(rep, n)),
              result)
  
  
  results_list[[row_index]] <- result
  
  # Close the database connection
  dbDisconnect(db)
}



## Check the result list & data structure
# length(results_list) 
# str(results_list[[1]])


## Combine the list of results and adjust the format of columns
df_LB <- bind_rows(results_list)

# unify time unit
df_LB$t <- trunc(df_LB$t) 

# convert p_profile from string to factor
df_LB$p_profile <- string_to_factor(df_LB$p_profile)

## Save combined data into a .csv file
data_path <- "/Users/YingJie/Documents/GitHub/ABM_code/analysis/df_LB_ctl_t20000.csv"
write.csv(df_LB, file = data_path, row.names = FALSE)