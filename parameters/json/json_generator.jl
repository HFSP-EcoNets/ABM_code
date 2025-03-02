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

# Generate same tensor as the previous function (the difference is after #/ Convert SparseArray to a format suitable for JSON)
function generate_propensities_json2(params, filename::AbstractString; power::Float64 = 2.0)
    #/ To compute the propensity tensor P, we need a list of all different strain ID & plasmid
    #  profiles that can exist in the system, even those that are not initially present.
    #!Note: The order of plasmid profiles is determined by Iterators.product.
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

end # module JSONgenerator
#/ End module