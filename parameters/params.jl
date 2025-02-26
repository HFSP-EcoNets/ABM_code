#/ Start module
module Params

#/ Packages
using Symbolics
using Parameters
using JSON

### FUNCTIONS
function load_parameters(filename::String)
	  file = JSON.read(filename, String)
    #/ Make symbolic dictionary from JSON file
    paramdict = Dict(
        (Symbol(key), value) for (key, value) in JSON.parse(file)
    )

    #/ As JSON does not contain strong types, we need to manually set the types here
    #@TODO Make this more robust/easier [low priority]
    competition_vectors = paramdict[:A] # this will give you Vectors of type Any
    for i in 1:length(competition_vectors) # convert every vector into type Float64
        competition_vectors[i] = Float64.(competition_vectors[i])
    end
    
    matrix_A = hcat(competition_vectors...) # combine vectors into a matrix
    paramdict[:A] = matrix_A # update the value of the symbol "A"

    compatibility_vectors = paramdict[:P]
    for i in 1:length(compatibility_vectors) # convert every vector into type Float64
        compatibility_vectors[i] = Float64.(compatibility_vectors[i])
    end

    infection_vectors = paramdict[:I] # this will give you Vectors of type Any
    for i in 1:length(infection_vectors) # convert every vector into type Float64
        infection_vectors[i] = Float64.(infection_vectors[i])
    end
    
    matrix_I = hcat(infection_vectors...) # combine vectors into a matrix
    paramdict[:I] = matrix_I # update the value of the symbol "I"
    
    matrix_P = hcat(compatibility_vectors...) # combine vectors into a matrix
    paramdict[:P] = matrix_P # update the value of the symbol "P"

    # keys "H" has values of type Any{Any}; let's convert this to Marix{Int64}
    HGT_vectors = paramdict[:H] # this will give you Vectors of type Any
    for i in 1:length(HGT_vectors) # convert every vector into type Float64
        HGT_vectors[i] = Int64.(HGT_vectors[i])
    end
    
    matrix_H = hcat(HGT_vectors...) # combine vectors into a matrix
    paramdict[:H] = matrix_H # update the value of the symbol "H"
    
    #/ Return parameter structure
    return params(; paramdict...)
end

### STRUCTURES
@with_kw mutable struct params
    "Simulation end time"
    t_final::Float64
    
    "Time between output events"
    t_output::Float64 = t_final
    
    "Seed"
    rng_seed::Union{Int64, Nothing}

    "Number of seeds (i.e. replicates)"
    n_seeds::Int

    "Number of initial bacterial strains = [n_b]"
    n_bstrains::Int

    "Number of initial plasmid strains = [n_p]"
    n_pstrains::Int
    
    "Number of bacterial substrains"
    n_bsubstrains::Int

    "Strain ID for each substrain"
    strain_id_for_each_substrain::Vector{Int}

    "Number of initial individuals of each bacterial substrain"
    n_ind_bsubstrains::Vector{Int}

    "Plasmid profiles of bacterial substrains"
    p_profile_bsubstrains::Vector{Vector{Int}}

    "Strain-specific maximum bacterial intrinsic growth rate (1/h) = [η_i]"
    growth_rate::Vector{Float64}

    "Constant strain bacterial carrying capacity (1/mL) = [K]"
    carrying_capacity::Vector{Float64}
    
    "Substrain-specific infection rate [γ_ij]"
    infection_rate::Vector{Float64}
    
    "Strain-specific probability of seggregation error [e_i]"
    segregation_error::Vector{Float64}
    
    "Strain-specific mininum bacterial intrinsic death rate (1/h) [μ_i]"
    death_rate::Vector{Float64}

    "Strain- and time- specific perturbation impact at t = 0, [ϵ(i)]"
    # unmutable; treated as amptitude when perturbaiton is fluctuating with time 
    perturbation_impact::Vector{Float64}

    "Strain-specific plasmid resistances = [r_β]"
    plasmid_resistance::Vector{Float64}

    "Strain-specific plasmid cost on bacterial growth = [c_β]"
    plasmid_cost::Vector{Float64}
    
    "Bacterial competition matrix = [A]"
    A::Matrix{Float64}

    "HGT matrix = [H]"
    H::Matrix{Int}
    
    "Infection matrix = [I]"
    I::Matrix{Float64}

    "Plasmid compatibility matrix = [P]"
    P::Matrix{Float64}
    
    "path to the propsensity JSON file"
    tensor_file::Union{AbstractString, Nothing}
end


end # module Params
#/ End module
