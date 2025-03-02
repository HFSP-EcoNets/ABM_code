### This is a Julia script that generate tensors locally for experimental design ###

#/ Laod Packages & module
using Revise
using CSV
using DataFrames
using MicrobePlasmidABM

#/ Assign the path to your CSV file (used only when running this script locally)
file_path = "input/3x4_test_IP_ctl_K_LB6-2_traits-fixed.csv"

#/ Reading the CSV file into a DataFrame and create a list of paths to the propensity JSON files (used only when running this script locally)
df = CSV.read(file_path, DataFrame)

#/ Loop through the list of paths to create the JSON files
#  Note: repliactes of the same experiements share the same tensor
rows = findall(x -> x == 1, df.rep)# rows that has experimental replicate = 1; assuming each experiment has a unique set of H, I, P matrices
for i in 1:rows
    path_parameter_file = df.ps[i]
    path_tensor_file = df.prop_tensor[i]
    params = Params.load_parameters(path_parameter_file);
    MicrobePlasmidABM.JSONgenerator.generate_propensities_json(params, path_tensor_file);                 
end    



