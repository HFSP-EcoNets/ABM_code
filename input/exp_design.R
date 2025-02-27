#### This is the R-script that generate a .csv table of YJ's test experimental design for a 3 x 4 system with all traits fixed across the strains & replicates, together with the .json files for simulations####
### compared to the exp_design_3x4_traits-fixed_YJ6.R, this R-script used t_final = 20000 rather than 2000.
### Download the packages (if packages are missing)
# install.packages("truncnorm")
# install.packages("jsonlite")
# install.packages("tidyverse")


### Load the packages & source script
library(truncnorm) # for sampling from a truncated normal distribution
library(jsonlite) # for converting a data frame to JSON format
library(tidyverse)
source("/Users/YingJie/Documents/GitHub/ABM_code/exp_design_util.R") #customized functions for experimental design

## The table contains multiple rows (i.e. simulations), with the following columns (*: vector of vector; stored as strings):
# 1.  job_id: the id for each simulation (i.e. job) ; an empty column that will only be filled while running the simulation pipeline
# 2.  key: the pre-assigned, never reused integer index of each simulation; for the job manager to locate the related .json files for simulations
# 3.  exp: the experiment ids for a set of simulations that has the same parameter set but different seeds for replicates
# 4.  ps: name of .json file that contains the related parameters
# 5.  label: the abstract profile of the simulation, e.g. A_std1_s123_YJ; could be used later for the name of .json input (key + label) and .sqlite output (job_id + key + label)
# 6.  rep: replicate number (currently 1:3) for each experiment id
# 7.  seed: the seed used for each simulation; each experiment shall have the same set of seeds
# 8.  n_b: number of initial bacterial strains in the system
# 9.  n_p: number of initial plasmid strains in the system
# 10. str_id: strain ids of the subpouplations 
# 11. N_subp: initial abundances (i.e. number of individuals) of the subpopulations in the system*
# 12. pp: initial plasmid profiles of the subpopulations in the system*
# 13. gr: intrinsic growth rates bacterial populations in the system*
# 14. K: carrying capacity of bacterial populations in the system*
# 15. inf_coef: intrinsic infection coefficient of bacterial populations in the system*
# 16. seg_err: probabilities of segregation error of bacterial populations in the system*
# 17. dr: intrinsic death rates of bacterial populations in the system*
# 18. ptb: perturbation coefficients of bacterial populations in the system*
# 19. pr: plasmid resistences of plasmid strains in the system*
# 20. pc: plasmid costs of plasmid strains in the system*
# 21. A: competition matrix* 
# 22. H: HGT matrix*
# 23. I: infection matrix*
# 24. P: plasmid compatibility matrix*
# 25. tensor_file: the file path to the JSON file of promensity tensor

### Step 0: check & set working directory
getwd()
setwd("/Users/YingJie/Documents/GitHub/ABM_code/input")

### Step 1: define constants & parameters (manual)

## Define the seeds and other constants
set.seed(123) # set the seeds for reproducibility (optional)

n_rep <- 300 # set the replicates (i.e. number of seeds required for each experimental treatment)
random_seeds <- sample(1:3000, n_rep) # generate three random integer seeds (note: the default replace = FALSE)


t_f <- 20000.0
t_o <- 5
n_s <- 1 # this is set to 1 for the HPC pipeline, for this argument won't be used there...

## Define the scale of the system (i.e. number of bacteria and plasmid strains)
n_bstrains <- 3
n_pstrains <- 4


####### Note: with all experiments, fix all plasmid and bacterial traits across strains

# matrix I

# full
vectors_I_tr1 <- lapply(1:n_pstrains, function(x) rep(1, n_bstrains))
str_vectors_I_tr1 <- paste(vectors_I_tr1, collapse = ",")

# nested
vectors_I_tr2 <- list(c(1, 1, 1), c(1, 1, 0), c(1, 1, 0), c(1, 0, 0))
str_vectors_I_tr2 <- paste(vectors_I_tr2, collapse = ",")

# modular
vectors_I_tr3 <- list(c(0, 0, 1), c(0, 1, 1), c(1, 1, 0), c(1, 0, 0))
str_vectors_I_tr3 <- paste(vectors_I_tr3, collapse = ",")

# matrix P (note that diagonal elements are always 0)

# full
vectors_P_tr1 <- lapply(1:n_pstrains, function(x) rep(1, n_pstrains))
for (i in 1:n_pstrains){
    vectors_P_tr1[[i]][i] <- 0
}
str_vectors_P_tr1 <- paste(vectors_P_tr1, collapse = ",")

# modular
vectors_P_tr2 <- list(c(0, 1, 0, 0), c(1, 0, 0, 0), c(0, 0, 0, 1), c(0, 0, 1, 0))
str_vectors_P_tr2 <- paste(vectors_P_tr2, collapse = ",")

# nested
vectors_P_tr3 <- list(c(0, 1, 1, 1), c(1, 0, 0, 0),  c(1, 0, 0, 0), c(1, 0, 0, 0))
str_vectors_P_tr3 <- paste(vectors_P_tr3, collapse = ",")



param1 <- list(vectors_I_tr1, vectors_I_tr2, vectors_I_tr3)
param2 <- list(vectors_P_tr1, vectors_P_tr2, vectors_P_tr3)

n_pcomb <- length(param1) * length(param2)
n_row <- n_pcomb * n_rep

## Define the initial number of agents (i.e. subpopulations) and related parameters

# Note: Link-balanced - design p_profiles & abundances based on the structure of I 
# each treatment I will have its own inital state, with non-zero I cells having a given abundance of 10
# all potential single-plasmid infected subpouplations are present and have the same abundance = 10, which is different from the plasmid-free subpopulation abundance = 2000

#total number of initial subpopulations (fixed across replicates)
num_subpops_total <- c(15,11,9)

#strain ids of initial subpopulations (fixed across replicates)
strain_ids_full <- c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3)
str_strain_ids_full <- paste(strain_ids_full, collapse = ",")

strain_ids_mod <- c(1,1,1,2,2,2,3,3,3)
str_strain_ids_mod <- paste(strain_ids_mod, collapse = ",")


strain_ids_nest <- c(1,1,1,1,1,2,2,2,2,3,3)
str_strain_ids_nest <- paste(strain_ids_nest, collapse = ",")

#plasmid profiles of initial subpopulations (designed accoridng to structure of I; fixed across replicates)

p_profiles_full <- list(c(0,0,0,0), c(1,0,0,0), c(0,1,0,0), c(0,0,1,0), c(0,0,0,1), c(0,0,0,0), c(1,0,0,0), c(0,1,0,0), c(0,0,1,0), c(0,0,0,1), c(0,0,0,0), c(1,0,0,0), c(0,1,0,0), c(0,0,1,0), c(0,0,0,1))
str_p_profiles_full <- paste(p_profiles_full, collapse = ",")

p_profiles_mod <- list(c(0,0,0,0), c(0,0,1,0), c(0,0,0,1), c(0,0,0,0), c(0,1,0,0), c(0,0,1,0), c(0,0,0,0), c(1,0,0,0), c(0,1,0,0))
str_p_profiles_mod <- paste(p_profiles_mod, collapse = ",")

p_profiles_nest <- list(c(0,0,0,0), c(1,0,0,0), c(0,1,0,0), c(0,0,1,0), c(0,0,0,1), c(0,0,0,0), c(1,0,0,0), c(0,1,0,0), c(0,0,1,0), c(0,0,0,0), c(1,0,0,0))
str_p_profiles_nest <- paste(p_profiles_nest, collapse = ",")

#initial subpopulation abundances (equal link abundace; abundance of infected individuals same among potential single-plasmid-infected subpopulations of each species (accroding to I design); fixed across repicates)
ini_abundances_full <- c(2000, 10, 10, 10, 10, 2000, 10, 10, 10, 10, 2000, 10, 10, 10, 10)
str_ini_abundances_full <- paste(ini_abundances_full, collapse = ",")

ini_abundances_mod <- c(2000, 10, 10, 2000, 10, 10, 2000, 10, 10)
str_ini_abundances_mod <- paste(ini_abundances_mod, collapse = ",")

ini_abundances_nest <- c(2000, 10, 10, 10, 10, 2000, 10, 10, 10, 2000, 10)
str_ini_abundances_nest <- paste(ini_abundances_nest, collapse = ",")

## Define the parameters for treatment (i.e. I and P) as string
I_levels_reps <- list()
P_levels_reps <- list()
for (i in seq_along(param1)){
    for (j in seq_along(param2)){
        I_level <- paste(param1[[i]], collapse = ",")
        P_level <- paste(param2[[j]], collapse = ",")
        for (j in 1:n_rep){
        I_levels_reps <- append(I_levels_reps, I_level)
        P_levels_reps <- append(P_levels_reps, P_level)
        }
    }
  
}

## Define the fixed strain-specific parameters, and matrices (as list of vector; each vector is a column of the matrices)

# growth rate
eta_i <- rep(1, n_bstrains) 
str_eta_i <- paste(eta_i, collapse = ",")

# death rate
mu_i <- rep(0.12, n_bstrains)
str_mu_i = paste(mu_i, collapse = ",")

# carrying capacity (fixed across strains & replicates)
K_i <- rep(2e4, n_bstrains) 
str_K_i = paste(K_i, collapse = ",")


# infection coefficient (fixed across strains & replicates)
gamma_i <- rep(1e-5, n_bstrains) 
str_gamma_i = paste(gamma_i, collapse = ",")


# probability of segregation error (fixed across strains & replicates)
e_i <- rep(1e-8, n_bstrains) 
str_e_i = paste(e_i, collapse = ",")


# perturbation coefficcient control (no perturbation; fixed across replicates)
epsilon_i_ctl <- rep(0.0, n_bstrains)
str_epsilon_i_ctl = paste(epsilon_i_ctl, collapse = ",")

# perturbation coefficient
epsilon_i <- rep(9, n_bstrains)
str_epsilon_i = paste(epsilon_i, collapse = ",")


# plasmid resistance
rho_i <- rep(0.7, n_pstrains)
str_rho_i = paste(rho_i, collapse = ",")


# plasmid cost
c_i <- rep(0.3, n_pstrains)
str_c_i = paste(c_i, collapse = ",")




# matrix A  (homogeneous; all off-diagonal elements = 0.01; all diagonal elements = 1.0):  
vectors_A <- lapply(1:n_bstrains, function(x) rep(0.01, n_bstrains))
for (i in 1:length(vectors_A)){
  vectors_A[[i]][i] = 1.0
}
# do.call(cbind, vectors_A) # convert vector of vectors into matrix

str_vectors_A = paste(vectors_A, collapse = ",")


# matrix H (homogeneous ; all elements = 1):  
vectors_H <- lapply(1:n_bstrains, function(x) rep(1, n_bstrains))
str_vectors_H = paste(vectors_H, collapse = ",")



### Step 2: Create the experimental design .csv  file in the target directory

## get the folder path where the experimental design .csv files are stored
folder_path <- "/Users/YingJie/Documents/GitHub/ABM_code/input"

## Create the data frame of experimental design (i.e. the .csv file)

# prepare the columns
lkey <- last_key(folder_path) # find the last key

job_id <- rep(NA, n_row)
key <- seq(lkey+1, lkey + n_row, by = 1)
exp <- rep(c(1:n_pcomb), each = n_rep)
ps <- c() # path to the parameter json files; to be filled using for-loop
prop_tensor <- c() # path to the propensity tensor json files; to be filled using for-loop
label <- c() # to be filled using for-loop; name of label should be e.g. gr_value_seed_Name
rep <- rep(1:n_rep, times = n_pcomb)
seed <- rep(random_seeds, times = n_pcomb)
n_b <- rep(n_bstrains, n_row)
n_p <- rep(n_pstrains, n_row)
#str_id <- rep(str_strain_ids, n_row)# for option 1 of initial state
str_id <- c(rep(str_strain_ids_full, n_row/3), rep(str_strain_ids_nest, n_row/3), rep(str_strain_ids_mod, n_row/3)) # for option 2 of initial state
#N_subp <- rep(str_ini_abundances, n_row) # for option 1 of initial state
N_subp <- c(rep(str_ini_abundances_full, n_row/3), rep(str_ini_abundances_nest, n_row/3), rep(str_ini_abundances_mod, n_row/3)) # for option 2 of initial state
#pp <- rep(str_p_profiles, n_row) # for option 1 of initial state
pp <- c(rep(str_p_profiles_full, n_row/3), rep(str_p_profiles_nest, n_row/3), rep(str_p_profiles_mod, n_row/3))  # for option 2 of initial state 
gr <- rep(str_eta_i, n_row) 
K <- rep(str_K_i, n_row)
inf_coef <- rep(str_gamma_i, n_row)
seg_err <- rep(str_e_i, n_row) 
dr <- rep(str_mu_i, n_row)
ptb <- rep(str_epsilon_i_ctl, n_row) # for groups without perturbation
#ptb <- rep(str_epsilon_i, n_row) # for gorups with perturbation
pr <- rep(str_rho_i, n_row)
pc <- rep(str_c_i, n_row)
A <- rep(str_vectors_A, n_row)
H <- rep(str_vectors_H, n_row)
I <- unlist(I_levels_reps) # the parameter for experimental treatment
P <- unlist(P_levels_reps) # the parameter for experimental treatment

I_labels <- c(rep("tr1", n_row/3), rep("tr2", n_row/3), rep("tr3", n_row/3))

P_labels <- c(rep("tr1", n_row/9), rep("tr2", n_row/9), rep("tr3", n_row/9),
              rep("tr1", n_row/9), rep("tr2", n_row/9), rep("tr3", n_row/9), 
              rep("tr1", n_row/9), rep("tr2", n_row/9), rep("tr3", n_row/9))


# set label, ps, and prop_tensor
# for groups without perturbation
for (i in 1:n_row){
  label_sub <- paste("ctl","I", I_labels[i], "P", P_labels[i],"seed", seed[i], "YJ", sep = "_")
  label <- append(label, label_sub)
  
  ps_sub <- paste("input/json/ctl_I", I_labels[i], "P", P_labels[i], "seed", seed[i], "YJ.json", sep = "_")
  ps <- append(ps, ps_sub)

  prop_tensor_sub <- paste("input/json/tensor_I", I_labels[i], "P", P_labels[i], "YJ.json", sep = "_")
  prop_tensor <- append(prop_tensor, prop_tensor_sub)
  
}

# for groups with perturbation
for (i in 1:n_row){
  label_sub <- paste("I", I_labels[i], "P", P_labels[i],"seed", seed[i], "YJ", sep = "_")
  label <- append(label, label_sub)
  
  ps_sub <- paste("input/json/I", I_labels[i], "P", P_labels[i], "seed", seed[i], "YJ.json", sep = "_")
  ps <- append(ps, ps_sub)

  prop_tensor_sub <- paste("input/json/tensor_I", I_labels[i], "P", P_labels[i], "YJ.json", sep = "_")
  prop_tensor <- append(prop_tensor, prop_tensor_sub)
  
}

# combine the columns to make the data frame
design <- tibble(
  job_id, key, exp, ps, prop_tensor, label, rep,
  seed, n_b, n_p, str_id, N_subp, pp, gr,
  K, inf_coef, seg_err, dr, ptb, pr,
  pc, A, H, I, P
)

# specify the file path where you want to save the CSV file

# for groups without perturabtion
file_path_name <- "/Users/YingJie/Documents/GitHub/ABM_code/input/3x4_test_IP_ctl_K_LB6-2_traits-fixed.csv"  # Lab's desktop's path; replace with the desired file path and name

# for groups with perturabtion
# file_path_name <- "/Users/YingJie/Documents/GitHub/ABM_code/input//3x4_test_IP_K_LB6-2_traits-fixed.csv"

# Write the data frame to a CSV file
write.csv(design, file_path_name, row.names = FALSE)




### Step 3: Create parameter JSON files in the target directory

## Set directory file path
directory_path <- "/Users/YingJie/Documents/GitHub/ABM_code/input/json"

## Create list, convert to and save the json files in the desired directory

# for Option 2 of initial state
for (i in 1:n_row){
  if (I_labels[i] == "tr1"){
    vectors_I <- vectors_I_tr1
    n_bsubstrains = num_subpops_total[1]
    strain_id_for_each_substrain = strain_ids_full
    n_ind_bsubstrains = ini_abundances_full
    p_profile_bsubstrains = p_profiles_full
    } 
  else if (I_labels[i] == "tr2") {
    vectors_I <- vectors_I_tr2
    n_bsubstrains = num_subpops_total[2]
    strain_id_for_each_substrain = strain_ids_nest
    n_ind_bsubstrains = ini_abundances_nest
    p_profile_bsubstrains = p_profiles_nest}
  else {
    vectors_I <- vectors_I_tr3
    n_bsubstrains = num_subpops_total[3]
    strain_id_for_each_substrain = strain_ids_mod
    n_ind_bsubstrains = ini_abundances_mod
    p_profile_bsubstrains = p_profiles_mod}

  if (P_labels[i] == "tr1"){
    vectors_P <- vectors_P_tr1} 
  else if (P_labels[i] == "tr2") {
    vectors_P <- vectors_P_tr2}
  else {
    vectors_P <- vectors_P_tr3}

  params <- list(
    t_final = t_f, 
    t_output = t_o,
    rng_seed = seed[i],
    n_seeds = n_s,
    n_bstrains = n_bstrains,
    n_pstrains = n_pstrains,
    n_bsubstrains = n_bsubstrains,
    strain_id_for_each_substrain = strain_id_for_each_substrain,
    n_ind_bsubstrains = n_ind_bsubstrains,
    p_profile_bsubstrains = p_profile_bsubstrains,
    growth_rate = eta_i,
    carrying_capacity = K_i,
    infection_rate = gamma_i,
    segregation_error = e_i,
    death_rate = mu_i,
    perturbation_impact = epsilon_i_ctl, # for groups without perturbation
    #perturbation_impact = epsilon_i, # for groups with perturbation
    plasmid_resistance = rho_i,
    plasmid_cost = c_i,
    A = vectors_A,
    H = vectors_H,
    I = vectors_I,
    P = vectors_P,
    tensor_file = prop_tensor[i]
  )
  json_output <- toJSON(params, pretty = TRUE, auto_unbox = TRUE)
  file_path <- file.path(directory_path, ps[i])
  writeLines(json_output, file_path)
  #cat("File", ps[i], "created.\n")
}
