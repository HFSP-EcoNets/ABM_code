### This is a Julia script that generate tensors locally for experimental design ###

#/ Laod Packages & module
using Revise
using CSV
using DataFrames
using MicrobePlasmidABM

#/ Assign the path to your CSV file (used only when running this script locally)
file_path = "hpc/sims/input/test_IP_design_I_ctl.csv"
file_path = "hpc/sims/input/test_IP_design_II_ctl.csv"
file_path = "hpc/sims/input/10x10_test_design.csv"

#/ Reading the CSV file into a DataFrame and create a list of paths to the propensity JSON files (used only when running this script locally)
df = CSV.read(file_path, DataFrame)

#/ Set the folder path to the json files
folder_path = "hpc/sims/"

#/ Loop through the list of paths to create the JSON files
#  Note: repliactes of the same experiements share the same tensor
rows = findall(x -> x == 1, df.rep)# rows that has experimental replicate = 1
for i in 1:rows
    path_parameter_file = df.ps[i]
    path_tensor_file = df.prop_tensor[i]
    params = Params.load_parameters(joinpath(folder_path, path_parameter_file));
    MicrobePlasmidABM.JSONgenerator.generate_propensities_json(params, joinpath(folder_path, path_tensor_file));                 
end    



