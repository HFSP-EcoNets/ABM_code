#= Create series of .json files for coputational experiment =#
#/ Start module
module JSONgenerator

#/ Packages
using JSON
using Base.Filesystem
using LinearAlgebra
using StaticArrays
using SparseArrayKit
using Dates

# Function for generating the JSON file of propensity tensor
function generate_propensities_json(params, filename::AbstractString, power::Float64 = 2.0)
    start_time = now()
    combs = Iterators.product([[0,1] for i in 1:params.n_pstrains]...) #/ create tuple (i.e. immutable collection of elements) of vectors[0,1] repeated by n_pstrains
    all_profiles = collect.(vec([comb for comb in combs])) #/ create all possible plasmid profiles in the system
    np = length(all_profiles) #/ number of all possible plasmid profiles in the system = 2^(n_pstrains)
    nb = params.n_bstrains #/ number of all possible strain profiles (i.e. strain IDs) in the system
    rep_strains = repeat(collect(1:nb), inner=np)
    rep_profiles = repeat(all_profiles, outer=nb)
    
    #/ Allocate propensity tensor 
    Γ = SparseArrayKit.SparseArray{Float64}(undef, np * nb, np * nb, np * nb)
    #/ Loop through all possible strain & plasmid profiles and infection pathways
    for (donor, (focal_strain, focal_plasmid_profile)) in enumerate(zip(rep_strains, rep_profiles)) #/ go through all potential donors' strain & plasmid profiles
        @info "loop at donor of strain & plasmid profile" donor focal_strain, focal_plasmid_profile 
        if any(x -> x > 0, focal_plasmid_profile) #/ if the donor has any plasmid 
            #@info "loop at donor of strain & plasmid profile" donor focal_strain, focal_plasmid_profile 
            for (recipient, (other_strain, other_plasmid_profile)) in enumerate(zip(rep_strains, rep_profiles)) # go through all potential recipients' strain & plasmid profiles
                #/ Check whether infection can occur from the HGT matrix H
                if  size(params.H) == (1,1) || params.H[other_strain, focal_strain] != 0
                    #/ Determine the no. of different plasmids that could be transferred
                    transferrable_plasmids = focal_plasmid_profile .> other_plasmid_profile
                    transferrable_plasmids = (transferrable_plasmids .* (params.I[other_strain,:] .> 0)) # filter out the plasmids that cannnot infect the recipient strain based on matrix I
                    #/ Get the indices (ids) of the plasmids that could be transferred
                    indices = findall(transferrable_plasmids)
                    #/ If there are transferrable plasmids (based on both the vacancy and infection matrix I), check all combinations
                    if any(transferrable_plasmids)
                        #/ Filter the transferr
                        #/ Gather all combinations of transferrable plasmids (i.e. possible outcomes of transfferable plasmid profiles of the transconjugant)
                        _combs = Iterators.product([[0,1] for p in 1:sum(transferrable_plasmids)]...)
                        _combs = collect.(vec(collect(_combs)))
                        #/ Loop through all combinations (expect for the first one, as it has none)
                        for comb in _combs[2:end]
                            #/ Gather infection profile
                            infection_profile = zeros(Int, params.n_pstrains) #/ create plasmid-free plasmid profile [0, 0, ...]
                            infection_profile[indices] .= comb #/ replace elements of plasmid profile with possible outcomes focusing on the transfferable plasmids 
                            #/ Compute the tranconjugant's profile after infection with infection profile
                            new_profile = other_plasmid_profile .+ infection_profile
                            #/ Check what profile that is (i.e. find the transconjugant's index)
                            transconjugant = findfirst(xy -> first(xy) == other_strain && last(xy) == new_profile, collect(zip(rep_strains, rep_profiles)))
                            
                            #/ Check element condition with matrices I & P and define element (tensor element will only be created if condition is met)
                            _all_plasmids = findall(new_profile .> 0)
            
                            if  size(params.I, 2) == 1 || all(params.P[other_plasmid,plasmid] != 0 for plasmid in _all_plasmids for other_plasmid in _all_plasmids[_all_plasmids .!= plasmid])
                                #!Note: We assume here that infecting 'p' plasmids scales as 1/2^(p-1)
                                Γ[transconjugant,donor,recipient] = 1.0 / power^(sum(infection_profile)-1)
                            end    
                            
                        end
                        #/ Normalize propensity tensor so that it can be used as a weight vector
                        #!Note: To avoid division by 0, we only do the normalization for sum(Γ[:, donor, recipient]) != 0
                        if sum(Γ[:, donor, recipient]) != 0
                        Γ[:,donor, recipient] ./= sum(Γ[:, donor, recipient])
                        end
                    end     
                end
            end
        end
    end
    #/ Convert SparseArray to a format suitable for JSON
    values = collect(nonzero_values(Γ))
    indices = Tuple.(collect(nonzero_keys(Γ)))
    nzs = [(c[1],c[2],c[3],γ) for (c,γ) in zip(indices, values)]
    Γdict = Dict(
        "dims" => size(Γ),
        "nonzeros" => nzs
    )

    #/ Write JSON to a file
    open(filename, "w") do io
    JSON.print(io, Γdict)
    end

    end_time = now()
    elapsed_seconds = Dates.value(end_time - start_time)/1000.0
    @info "time elasped for this tensor generation (s):" elapsed_seconds 

    return Γ # optional
end

# previous function (the difference is after #/ Convert SparseArray to a format suitable for JSON)
function generate_propensities_json2(params, filename::AbstractString; power::Float64 = 2.0)
    #/ To compute the propensity tensor P, we need a list of all different strain ID & plasmid
    #  profiles that can exist in the system, even those that are not initially present.
    #!Note: The order of plasmid profiles is determined by Iterators.product.
    #@TODO: This ordering might prove problematic as, in principle, the input parameters
    #       can contain an order that is *not* realized by Iterators.product.
    #@TODO: Right now we can only define the propensity tensor here, but we would like it
    #       to be an input parameter to the system instead.
    combs = Iterators.product([[0,1] for i in 1:params.n_pstrains]...) #/ create tuple (i.e. immutable collection of elements) of vectors[0,1] repeated by n_pstrains
    all_profiles = collect.(vec([comb for comb in combs])) #/ create all possible plasmid profiles in the system
    np = length(all_profiles) #/ number of all possible plasmid profiles in the system = 2^(n_pstrains)
    nb = params.n_bstrains #/ number of all possible strain profiles (i.e. strain IDs) in the system
    rep_strains = repeat(collect(1:nb), inner=np)
    rep_profiles = repeat(all_profiles, outer=nb)
    
    #/ Allocate propensity tensor 
    Γ = SparseArrayKit.SparseArray{Float64}(undef, np * nb, np * nb, np * nb)
    #/ Loop through all possible strain & plasmid profiles and infection pathways
    for (donor, (focal_strain, focal_plasmid_profile)) in enumerate(zip(rep_strains, rep_profiles)) #/ go through all potential donors' strain & plasmid profiles
        @info "loop at donor of strain & plasmid profile" donor focal_strain, focal_plasmid_profile 
        if any(x -> x > 0, focal_plasmid_profile) #/ if the donor has any plasmid 
            #@info "loop at donor of strain & plasmid profile" donor focal_strain, focal_plasmid_profile 
            for (recipient, (other_strain, other_plasmid_profile)) in enumerate(zip(rep_strains, rep_profiles)) # go through all potential recipients' strain & plasmid profiles
                #/ Check whether infection can occur from the HGT matrix H
                if  size(params.H) == (1,1) || params.H[other_strain, focal_strain] != 0
                    #/ Determine the no. of different plasmids that could be transferred
                    transferrable_plasmids = focal_plasmid_profile .> other_plasmid_profile
                    transferrable_plasmids = (transferrable_plasmids .* (params.I[other_strain,:] .> 0)) # filter out the plasmids that cannnot infect the recipient strain based on matrix I
                    #/ Get the indices (ids) of the plasmids that could be transferred
                    indices = findall(transferrable_plasmids)
                    #/ If there are transferrable plasmids (based on both the vacancy and infection matrix I), check all combinations
                    if any(transferrable_plasmids)
                        #/ Filter the transferr
                        #/ Gather all combinations of transferrable plasmids (i.e. possible outcomes of transfferable plasmid profiles of the transconjugant)
                        _combs = Iterators.product([[0,1] for p in 1:sum(transferrable_plasmids)]...)
                        _combs = collect.(vec(collect(_combs)))
                        #/ Loop through all combinations (expect for the first one, as it has none)
                        for comb in _combs[2:end]
                            #/ Gather infection profile
                            infection_profile = zeros(Int, params.n_pstrains) #/ create plasmid-free plasmid profile [0, 0, ...]
                            infection_profile[indices] .= comb #/ replace elements of plasmid profile with possible outcomes focusing on the transfferable plasmids 
                            #/ Compute the tranconjugant's profile after infection with infection profile
                            new_profile = other_plasmid_profile .+ infection_profile
                            #/ Check what profile that is (i.e. find the transconjugant's index)
                            transconjugant = findfirst(xy -> first(xy) == other_strain && last(xy) == new_profile, collect(zip(rep_strains, rep_profiles)))
                            
                            #/ Check element condition with matrices I & P and define element (tensor element will only be created if condition is met)
                            _all_plasmids = findall(new_profile .> 0)
            
                            if  size(params.I, 2) == 1 || all(params.P[other_plasmid,plasmid] != 0 for plasmid in _all_plasmids for other_plasmid in _all_plasmids[_all_plasmids .!= plasmid])
                                #!Note: We assume here that infecting 'p' plasmids scales as 1/2^(p-1)
                                Γ[transconjugant,donor,recipient] = 1.0 / power^(sum(infection_profile)-1)
                            end    
                            
                        end
                        #/ Normalize propensity tensor so that it can be used as a weight vector
                        #!Note: To avoid division by 0, we only do the normalization for sum(Γ[:, donor, recipient]) != 0
                        if sum(Γ[:, donor, recipient]) != 0
                        Γ[:,donor, recipient] ./= sum(Γ[:, donor, recipient])
                        end
                    end     
                end
            end
        end
    end
    #/ Convert SparseArray to a format suitable for JSON
    # Γ_json = Dict("dims" => size(Γ), "nonzeros" => [(i, j, k, Γ[i, j, k]) for (i, j, k, _) in findall(!iszero, Γ)])
    Γ_json = Dict("dims" => size(Γ), "nonzeros" => [(ci[1], ci[2], ci[3], Γ[ci]) for ci in findall(!iszero, Γ)])
    #/ Write JSON to a file
    open(filename, "w") do io
    JSON.print(io, Γ_json)
    end

    return Γ # optional
end

function generate_propensities_json_simple(params, filename::AbstractString)
    #/ To compute the propensity tensor P, we need a list of all different strain ID & plasmid
    start_time = now()
    combs = Iterators.product([[0,1] for i in 1:params.n_pstrains]...) #/ create tuple (i.e. immutable collection of elements) of vectors[0,1] repeated by n_pstrains
    all_profiles = collect.(vec([comb for comb in combs])) #/ create all possible plasmid profiles in the system
    np = length(all_profiles) #/ number of all possible plasmid profiles in the system = 2^(n_pstrains)
    nb = params.n_bstrains #/ number of all possible strain profiles (i.e. strain IDs) in the system
    rep_strains = repeat(collect(1:nb), inner=np)
    rep_profiles = repeat(all_profiles, outer=nb)
    
    #/ Allocate propensity tensor 
    Γ = SparseArrayKit.SparseArray{Float64}(undef, np * nb, np * nb, np * nb)
    #/ Loop through all possible strain & plasmid profiles and infection pathways
    for (donor, (focal_strain, focal_plasmid_profile)) in enumerate(zip(rep_strains, rep_profiles)) #/ go through all potential donors' strain & plasmid profiles
        @info "loop at donor of strain & plasmid profile" donor focal_strain, focal_plasmid_profile 
        if any(x -> x > 0, focal_plasmid_profile) #/ if the donor has any plasmid 
            #@info "loop at donor of strain & plasmid profile" donor focal_strain, focal_plasmid_profile 
            for (recipient, (other_strain, other_plasmid_profile)) in enumerate(zip(rep_strains, rep_profiles)) # go through all potential recipients' strain & plasmid profiles
                #/ Check whether infection can occur from the HGT matrix H
                if  size(params.H) == (1,1) || params.H[other_strain, focal_strain] != 0
                    #/ Determine the no. of different plasmids that could be transferred
                    transferrable_plasmids = focal_plasmid_profile .> other_plasmid_profile
                    transferrable_plasmids = (transferrable_plasmids .* (params.I[other_strain,:] .> 0)) # filter out the plasmids that cannnot infect the recipient strain based on matrix I
                    #/ Get the indices (ids) of the plasmids that could be transferred
                    indices = findall(transferrable_plasmids)

                    # If recipient is not plasmid-free, filter out the plasmids that cannnot infect the recipient based on matrix P
                    if sum(other_plasmid_profile) != 0
                        other_plasmids = findall(x -> x != 0, other_plasmid_profile)
                        for plasmid in indices
                           if any(params.P[other_plasmid,plasmid] == 0 for other_plasmid in other_plasmids[other_plasmids .!= plasmid])
                            transferrable_plasmids[plasmid] = false 
                           end 
                        end  
                    end    
                      
                    #/ If there are transferrable plasmids (based on vacancy, infection matrix I, and plasmid compatibiltiy matrix P), check all combinations
                    if any(transferrable_plasmids)
                        #/ Update the plasmid indice of transferrable_plasmids
                        indices = findall(transferrable_plasmids)
                        #/ Gather all combinations of transferrable plasmids (i.e. possible outcomes of transfferable plasmid profiles of the transconjugant)
                        _combs = Iterators.product([[0,1] for p in 1:sum(transferrable_plasmids)]...)
                        _combs = collect.(vec(collect(_combs)))
                        #/ Exclude combinations that has no or more than one plasmid
                        _combs = [vec for vec in _combs if sum(vec) == 1]
                        #/ Loop through all combinations 
                        for comb in _combs
                            #/ Gather infection profile
                            infection_profile = zeros(Int, params.n_pstrains) #/ create plasmid-free plasmid profile [0, 0, ...]
                            infection_profile[indices] .= comb #/ replace elements of plasmid profile with possible outcomes focusing on the transfferable plasmids 
                            #/ Compute the tranconjugant's profile after infection with infection profile
                            new_profile = other_plasmid_profile .+ infection_profile
                            #/ Check what profile that is (i.e. find the transconjugant's index)
                            transconjugant = findfirst(xy -> first(xy) == other_strain && last(xy) == new_profile, collect(zip(rep_strains, rep_profiles)))
                            #/ Fill in the non-zero tensor element 
                            Γ[transconjugant,donor,recipient] = 1.0 
                            
                        end
                        #/ Normalize propensity tensor so that it can be used as a weight vector
                        #!Note: To avoid division by 0, we only do the normalization for sum(Γ[:, donor, recipient]) != 0
                        if sum(Γ[:, donor, recipient]) != 0
                        Γ[:,donor, recipient] ./= sum(Γ[:, donor, recipient])
                        end
                    end     
                end
            end
        end
    end

    #/ Convert SparseArray to a format suitable for JSON
    values = collect(nonzero_values(Γ))
    indices = Tuple.(collect(nonzero_keys(Γ)))
    nzs = [(c[1],c[2],c[3],γ) for (c,γ) in zip(indices, values)]
    Γdict = Dict(
        "dims" => size(Γ),
        "nonzeros" => nzs
    )

    #/ Write JSON to a file
    open(filename, "w") do io
    JSON.print(io, Γdict)
    end

    end_time = now()
    elapsed_seconds = Dates.value(end_time - start_time)/1000.0
    @info "time elasped for this tensor generation (s):" elapsed_seconds

    return Γ
end

# Function for generating the propensity tensor from its JSON file
function recreate_propensities_from_json(filename::AbstractString)
    # Read and parse the JSON file
    raw_json = open(filename, "r") do file
        read(file, String)
    end
    json_data = JSON.parse(raw_json)

    # Extract dimensions and non-zero elements
    dims = tuple(json_data["dims"]...)
    nonzeros = json_data["nonzeros"]

    # Initialize a sparse array with the given dimensions
    Γ_recreated = SparseArrayKit.SparseArray{Float64, 3}(undef, dims)

    # Populate the sparse array with non-zero elements
    for elem in nonzeros
        i, j, k, value = elem
        Γ_recreated[i, j, k] = value
    end

    return Γ_recreated
end

# Function for parameter sweep: e.g. sweeping over infection rate and plasmid cost
function jsons(params_1, params_2)
    #/ Define the range of parameters for parameter sweep
    #plasmid_cost_values = 0.0:0.05:1.0
    #encounter_rate_values = LinRange(1e-6, 1e-4, length(plasmid_cost_values))


    #/ Directory where JSON files will be saved
     #output_directory = "/Users/YingJie/Documents/GitHub/population_dynamics_ABM/src/parameters/json/sweep"
     output_directory = "/Users/lanexran/Documents/GitHub/population_dynamics_ABM/src/parameters/json/sweep" 


    #/ Create the output directory if it doesn't exist
    if !isdir(output_directory)
        mkdir(output_directory)
        #mkpath(output_directory)
    end


    #/ Loop through parameter values and create JSON files
    for i in 1:length(params_1)
        for j in 1:length(params_2)
        # Define the content of the JSON object with the current parameters
            json_content = Dict(
            "t_final" => 500.0,
            "t_output" => 5.0,
            "rng_seed" => nothing,
            "n_seeds" => 1,
            "n_bstrains" => 2,
            "n_pstrains" => 1,
            "n_bsubstrains" => 2,
            "strain_id_for_each_substrain" => [1, 1],
            "n_ind_bsubstrains" => [1e3, 1e2],
            "p_profile_bsubstrains" => [[0], [1]],
            "growth_rate" => [1.0],
            "death_rate" => [0.05],
            "carrying_capacity" => [1e4],
            "infection_rate" => [params_2[j]],  # Set the "encounter_rate" parameter
            "segregation_error" => [0.0],
            "perturbation_impact" => [0.0],
            "plasmid_resistance" => [1.0],
            "plasmid_cost" => [params_1[i]],  # Set the "plasmid_cost" parameter
            "A" => [[1.0]],
            "H" => [[1]],
            "I" => [[1.0]],
            "P" => [[0.0]],
            "tensor_file": nothing
            )

            # Generate the JSON file name based on the parameters
            json_filename = joinpath(output_directory, "param1_$(i)_param2_$(j).json")

            # Serialize the JSON content to a file
            open(json_filename, "w") do file
            JSON.print(file, json_content)
            end
        end
    end
end


# Function for testing distributions of matrix A (i.e. Var[A_ij]) with n bacterial strains and 0 plasmid
function jsons_A(nbstrains, nsample, std::Float64, mean::Float64=0.0, c::Float64=0.5)
    # nsample: statistically efficient sample size = 30
    # mean[A_ij] = 0
    # std = standard deviation = σ = 1.0 vs. 0.5
    # Var[A_ij] = cσ^2/nbstrains   

    #/ Directory where JSON files will be saved
     output_directory = "/Users/YingJie/Documents/GitHub/population_dynamics_ABM/src/parameters/json/sweep" # for lab laptop
     #output_directory = "/Users/lanexran/Documents/GitHub/population_dynamics_ABM/src/parameters/json/sweep" # for personal laptop

    #/ Create the output directory if it doesn't exist
    if !isdir(output_directory)
        mkdir(output_directory)
        #mkpath(output_directory)
    end

    #/create vectors
    strain_ids = collect(1:nbstrains)
    abundances = fill(1000, nbstrains)
    p_profiles = [fill(0, 1) for _ in 1:nbstrains]
    growth_rate = fill(1, nbstrains)
    death_rate = fill(0.05, nbstrains)
    carrying_capacity = fill(1e4, nbstrains)
    infection_rate = fill(1e-5, nbstrains)
    segregation_error = fill(0, nbstrains)
    perturbation_impact = fill(0.0, nbstrains)

    H = [fill(1, nbstrains) for _ in 1:nbstrains]
    I = fill(1, nbstrains)

    

    #/ Loop through parameter values and create JSON files
    for i in 1:nsample
        # Create matrix A with disried mean and distribution + based on radom structure of connectance (negatively correlated with phylogenetic distance), with A_ii = 1 and A_ij = 0 if connectance < 0.5
        # (Note: following Han's method) 
        a = abs.((std/sqrt(nbstrains)) .* randn(nbstrains,nbstrains) .+ mean/nbstrains) # generate a rondom number from a standard normal distribution (mean 0 and std = 1), multiply elements by std/sqr(nbstrains), add by mean/nbstrains, and take absolute values
        conn = rand(nbstrains, nbstrains)
        tconn = tril(conn)
        tconn[diagind(tconn)] .= 0.0
        conn = tconn + tconn'
        a[conn .< c] .= 0.0
        a[diagind(a)] .= 1

        # Convert A from matrix into vector of vectors (each vector is a column of A) for the JSON matrix format
        num_columns = size(a,2)
        A = [a[:, j] for j in 1:num_columns]

        # Define the content of the JSON object with the current parameters
        json_content = Dict(
        "t_final" => 500.0,
        "t_output" => 5.0,
        "rng_seed" => nothing,
        "n_seeds" => 1,
        "n_bstrains" => nbstrains,
        "n_pstrains" => 1,
        "n_bsubstrains" => nbstrains,
        "strain_id_for_each_substrain" => strain_ids,
        "n_ind_bsubstrains" => abundances,
        "p_profile_bsubstrains" => p_profiles,
        "growth_rate" => growth_rate,
        "death_rate" => death_rate,
        "carrying_capacity" => carrying_capacity,
        "infection_rate" => infection_rate, 
        "segregation_error" => segregation_error,
        "perturbation_impact" => perturbation_impact,
        "plasmid_resistance" => [0.0],
        "plasmid_cost" => [0.0],
        "A" => A,
        "H" => H,
        "I" => I,
        "P" => [[0.0]]
        )

        # Generate the JSON file name based on the parameters
        json_filename = joinpath(output_directory, "random_A_std$(std)_$(i).json")

        # Serialize the JSON content to a file
        open(json_filename, "w") do file
        JSON.print(file, json_content)
        end
        
    end
end

end # module JSONgenerator
#/ End module