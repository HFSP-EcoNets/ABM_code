# ABM code
This repository contains the ABM code for publication.

## :open_file_folder: Directory overview
Below an overview of repositories, codes and files used for simulation & data processing practices; only the ones that are currently used in the latest ABM are listed (the rest repositories and files are to be re-constructed/deprecated):

```bash
.
├── README.md
├── data
├── analysis
│   ├── Appendix_B.Rmd
│   ├── figure_themes.R
│   ├── process_util.R
│   └── gatheroutput.jl
├── gillespie
│    ├── gillespie.jl
│    ├── gillespie1.jl
│    ├── system.jl
│    └── system1.jl
├── parameters
│   ├── json
│   │   ├── json_generator.jl
│   │   ├── parameters_example3.json
│   │   ├── tensor_example3.json
│   │   └── tensor_json_generator.jl
│   └── params.jl 
├── utils
│   └── utils.jl
├── input
│   │   ├── 3x4_test_IP_ctl_K_LB6-2_traits-fixed.csv
│   │   ├── exp_design_util.R
│   │   └── exp_design.R
│   └── json
├── output
└── MicrobePlasmidABM.jl
```

## :warning: Advised workflow 
* Create a new branch for yourself with some descriptive new branch name.
* Make sure you're under your own branch before working on the codebase.
* Always pull from your branch on the cloud (i.e. source) before working on the local scripts/files.
* When satisfied, save, commit, and push the local adjustments onto the cloud.
* If needed, you may rebase your branch, or merge another branch into your branch. 
* Consult with senior reserachers before pulling/pushing/rebasing/merging if you're not sure how to do it.



## :cd: Installation (under Visual Studio Code)
### :package: Importing and installing the Julia project 
If you haven't included the local module (i.e. package) `MicrobePlasmidABM` in your VSCode environment, activate and instantiate the module following the below instructions:

1. Clone the repository `ABM_code` 
2. Open the Julia REPL terminal (by pressing "option" + "o" + "j" in the TERMINAL panel)
3. Under the Julia REPL terminal (julia>), run `cd("/path/to/ABM_code")`
4. Press "]" to change into pkg>, run `activate .` then run `instantiate`
5. Press "backspace" to return from pkg> to julia>, run  `using Pkg` then `Pkg.instantiate()` 
6. Run `using MicrobePlasmidABM` to load the module
7. Advice: run `import Pkg; Pkg.add("Revise")` and `using Revise`. This package traces and recompiles the changes of the module files within the REPL, so one doesn't need to reload Julia REPL after revising the code.
8. Optional: run `Pkg.resolve()` to determine the complete set of packages and their versions that satisfy the constraints specified in the project environment; run `Pkg.update()` update the global packages.

All the required packages should be downloaded and installed. 
The initial load might take some time, but subsequent loads should be faster.

After the installation is completed, simply run `using Revise` and `using MicrobePlasmidABM` next time in Julia REPL terminal for simulations.



## :books: Introduction of the Microbe-plasmid agent-based model (ABM)
For detailed information about the mathematical formulation, please read the supplemetary materials Appendix A of the manuscript:

(To do: add the link to the preprint/publication here)

### Model structure and simulation algorithm

This ABM treats each bacterial subpopulation (i.e. substrain) as an agent (i.e. an entity). A subpopulation is characterized by its strain i and plasmid plasmid profile (i.e. the plasmids it hosts) p.

With each bacterial strain, each subpopulation has a plasmid proflie, a vector of length = number of plasmid strains in the system. For example, in a community with two plasmid strains, the plasmid profile of a bacterial subpopulation will be vector of length = 2. The element(s) of the plasmid profile represent the presence (0) or absence (1) of the plasmid strain. Therefore, a system with one bacterial strain and n_p plasmid strains can have at most 2^n infection states (maximum length (p) = 2^n_p) and therefore 2^n_p subpopulations. For a system with one bacterial strain and two plasmid strains, the potential subpopulaions are listed as below.

i   p
1 [0,0],
1 [0,1],
1 [1,0],
1 [1,1]

Based on the strain, plasmid profile, and abundance, each subpopulation has its corresponding rates of actions (i.e. events) in the ABM.

The ABM is a stochastic, continuous-time, discrete-event simulation (Gillespie-style) that steps from event to event, with a random amount of time, real-valued, between events.

### Model state 

Model state is represented by the struct `states` (see `system.jl`), which contains the following fields:
* N: Total cell abundance across subpopulations
* abundances of subpopulations
* strain IDs of subpopulations  
* all potential plasmid profiles
* plasmid_profiles of subpouplations
* donors indices in the system
* donors indices in the propensity tensor
* recipients indices in the system (donor specific)
* recipient indices in the propensity tensor
* strain IDs in propensity tensor
* plasmid profiles in propensity tensor

Conceptually, it represents a collection of bacterial subpopulations and the plasmid strain(s) they host. The initial state of the system will be used to compute the rates of events, which will be further sued to sample the next time t and event that take place at that t.

### Model events and corresponding stage changes

There are five events for the bacterial subpopulations:
1. Growth: chosen subpopulation increases its abundance by 1.
2. Death: chosen subpopulatoin decreases its abundance by 1.
3. Segregatoin: the plasmid-free subpouplation with the same strain ID as the chosen subpopulation increases its abundance by 1.
4. Competition: chosen subpopulation decreases its abundance by 1.
5. Infection: chosen recipient subpopulation decreases its abundance by 1, while chosen transconjugant subpopulation increases its abundacne by 1.

For details about the equations and implementation of these rates, check the Overleaf manucript HFSP_ABM and the script `system.jl`.

### Model simulation algorithm
This ABM uses the Gillespie algorithm for simuation, with the follwing steps:

1. The perturbation impact is extracted based on the previous time step (for the first time step, this will be perturbation impact at time 0) 
2. The time to the next event is computed by drawing a waiting time from a Poisson process (i.e., exponentially distributed inter-event times). 
3. The type of event is chosen with weight proportional to its rate relative to the total rate of all events.
4. Under the chosen event, the specific subpopulation(s) is/are chosen sampled with weight(s) proportional to its abundance * event rate relative to the total abundance * event rate.
5. The corresponding state change(s) is/are executed accorind to the chosen event and subpopulation(s).
6. When simulation time meets the time to be recorded (e.g. t = 5, 10, 15, 20...), the time and corresponding state are temporarily recorded into some dataframes.
7. When simulation time meets the desired time steps for exporting data, the recorded data will be exported into a .sqlite ouput for futher processing using the R pipeline.

For details about the implementation of simulation algorithm, check the Overleaf manucript HFSP_ABM and the script `gillespie.jl`.

## :computer: Running simulations 

### Single simulation
You can run a single simulation given the location of the parameters file in JSON format, the parameters, and a single-line argument for running the simulatoin, such as:

```sh
filename = "src/parameters/json/parameters_example3.json"
params = Params.load_parameters(filename); 
output_folder::String = "data" 
job_key::String = "000"
JOB_ID::String = "000000"   
MicrobePlasmidABM.Gillespie.run(params; filename = filename, output_folder, job_key, JOB_ID)
```

After which, you shall see an automatically generated file `output000.sqlite` in the repository `data`.

### Mutiple simulations: replicates (sequential simulation on local computer; not needed with the HPC)
You can run mutiple simulations sequentially to generate a series of .sqlite output given the location of the parameters files in JSON format, the parameters, and a for-loop for running the simulatoins, such as:

```sh
job_keys = 17929:17931
job_seeds = [2463,2511,2227]
job_length = length(job_keys)   
output_folder::String = "data"
for i in 1:job_length
    job_key::String = string.(job_keys[i])
    JOB_ID::String = lpad(string(i), 6, "0")
    filename = "input/json/ctl_I_tr1_P_tr1_seed$(job_seeds[i])_YJ.json"
    params = Params.load_parameters(filename);
    MicrobePlasmidABM.Gillespie.run(params; filename = filename, output_folder, job_key, JOB_ID)
end 
```

After which, you shall see an automatically generated files `output$job_key.sqlite` in the repository `data`.

### Infection tensor generatoin
Each parameter input file (.json file) will contain a path to a infection tensor file (.json file) used for the simulation. This 3-dimentional tensor contains propensities at which a transconjugant of a certain subpopulation can be produed during an infection event, given the H, I, P networks, and the identities of plasmid donor and recipient. When not specifically asigned (i.e. "tensor_file": null), this tensor will be generated during the simulation. For computational efficiency, it is recommanded to generate the infeciton tensor file in advance using the parameter input file. For example, delete the tensor file `tensor_exmaple3.sqlite` in the repository `parameters/json`, then follow the instructions below:

```sh
# Set the folder path to the parameter input file
folder_path = "parameters/json/"

# Specificy the JSON file of parameters
filename = joinpath(folder_path, "parameters_example3.json")

# Get the parameters
params = Params.load_parameters(filename);

# Specificy the path and name of the JSON file of infection tensor
propensity_filename = params.tensor_file

# Generate the tensor file         
MicrobePlasmidABM.JSONgenerator.generate_propensities_json(params, propensity_filename); 
```

After which, you shall see an automatically generated files `tensor_exmaple3.sqlite` in the repository `parameters/json`.


### Input preparation and output analysis for mutiple simulations
Input (.json files) generation and output (.sqlite files) analysis is performed in R. For the R script of input generation, please check `exp_design.R` in the repository `input`. For the R script of output analysis, please check `list_output.R` and `Appendix_B.Rmd` in the repository `analysis`.

(To do: upload input & output of emperically based simulations and Appendix_C.Rmd if needed.)

### Input and output used in the manuscript
The raw input (.json files) and output (.sqlite files) used in the manuscript are respectively stored in the repository `input/json` and repository `output`. In addition, an overview of the input is also provided in .csv files, respectively stored in the repository `input`. 

## :page_facing_up: Notes 
1. Setting a high community-wise carrying capacity (e.g. 1e10) may greatly slows down the simulation as delta_t becomes extremely small in each loop. Therefore it is recommended to scale K below 1e6 for computational tractability.

