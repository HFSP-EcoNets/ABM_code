#= Specify the agent-based system using in the Gillespie algorithm =#
#= Pruturbatoin mode: single pulse, constant=#
#/ Start module
module System1

#/ Packages
using StaticArrays
using SparseArrayKit
using Dates
using LinearAlgebra

#/ Local packages
using MicrobePlasmidABM

### FUNCTIONS
function initialize_system(time, params)
    #/ generate propensitity tensor (i.e. Γ)
    @info "Getting propensities..."
    propensities = (params.tensor_file === nothing) ?
    get_propensities3(params) : MicrobePlasmidABM.JSONgenerator.recreate_propensities_from_json(params.tensor_file); 

    #/ Initialize the substrain arrays
    @info "Initializing strains..."
	states = initialize_strains(propensities, params)


    @info "Getting rates..."
    rates = get_rates(time, states, params)
    
    #/ Extract rates
    return states, rates, propensities
end

### STRUCTURES
# We would like the states and rates to be packed into a structure that is easy
# to pass around, instead of having functions with many arguments. The main reason
# for this is because many elementary reactions can change the state arrays of the
# system (e.g. when a subpopulation goes extinct, it needs to be removed from all
# relevant arrays). Doing this with seperate arrays is quite cumbersome, so instead
# we use a struct.
mutable struct states
    N::Int                                     # Total number of subpopulations
    abundances::Array{Int, 1}                  # Array of integer abundances
    subpopulations::Array{Int, 1}              # Array of strain indices
    all_profiles::Array{Array{Int, 1}, 1}      # Array of array of all plasmid profiles
    plasmid_profiles::Array{Array{Int, 1}, 1}  # Array of array of existing plasmid profiles
    donors::Array{Int, 1}                      # Array of list indices (in system) of donors
    donors_indices_Γ::Array{Int, 1}          # Array of list indices (in Γ) of donors
    recipients::Array{Array{Int, 1}, 1}        # Array of array of recipients (in system) for each donor
    recipients_indices_Γ::Array{Array{Int, 1}, 1} # Array of array of recipients (in Γ) for each donor
    rep_strains::Array{Int64, 1}               # Array of strain in Γ (fixed)
    rep_profiles::Array{Array{Int64, 1}, 1}    # Array of plasmid_profile in Γ (fixed)
end

mutable struct rates
    growth_rates::Array{Float64, 1}
    death_rates::Array{Float64, 1}
    segregation_rates::Array{Float64, 1}
    competition_rates::Array{Array{Float64, 1}, 1}
    infection_rates::Array{Float64, 1} # for the new infection equation
end

### HELPER FUNCTIONS
function initialize_strains(propensities, params)
    #/ Allocate arrays
    subpopulations = [strain_id for strain_id in params.strain_id_for_each_substrain] # strain id of subpopulations
    N = length(subpopulations) # number of subpopulations
    plasmid_profiles = [plasmid_profile for plasmid_profile in params.p_profile_bsubstrains]
    abundances = copy(params.n_ind_bsubstrains) # abundance of subpopulations
    #/ Get all possible plasmid profiles with a fixed ordering
    combs = Iterators.product([[0,1] for i in 1:params.n_pstrains]...)
    #~ Extract profiles by collecting the iterator and turning them into vectors
    all_profiles = collect.(vec(collect(combs)))
    #~ Collect the strains & plasmid profiles of subpulation in the propensity tensor that determined compartment flows.
    np = length(all_profiles) #/ number of all possible plasmid profiles in the system = 2^(n_pstrains)
    nb = params.n_bstrains #/ number of all possible strain profiles (i.e. strain IDs) in the system
    rep_strains = repeat(collect(1:nb), inner=np)
    rep_profiles = repeat(all_profiles, outer=nb)
    
    #/ Construct array of *indices* of the donors in the substrain array such that
    #  substrains[donors] will always return all possible donors in the system
    #!Note: Those with a plasmid are always donors
    donors = Array{Int}(undef, 0)
    donors_indices_Γ =  Array{Int}(undef, 0)
    for (i, (strain, plasmid_profile)) in enumerate(zip(subpopulations, plasmid_profiles))
        if any(plasmid_profile .> 0) 
            push!(donors, i)
            donor_index_Γ = findfirst(xy -> first(xy) == strain && last(xy) == plasmid_profile, collect(zip(rep_strains, rep_profiles)))
            push!(donors_indices_Γ, donor_index_Γ)
        end    
    end

    #/ Construct array of *indices* of the recipients for each donor such that,
    #  for each donor, we can obtain the recipients as
    #  recipients[donor]

    # make a list of valid recipients based on the propensity tensor
    recipients = Array{Array{Int, 1}, 1}(undef, 0)
    recipients_indices_Γ = Array{Array{Int, 1}, 1}(undef, 0)
    for i in donors
        #/ Get plasmid profile of donor
        donor_strain = subpopulations[i]
        donor_plasmid_profile = plasmid_profiles[i]
        recipients_of_donor = Array{Int, 1}(undef, 0)
        recipient_indices_Γ = Array{Int, 1}(undef, 0)

        donor_index_Γ = findfirst(xy -> first(xy) == donor_strain && last(xy) == donor_plasmid_profile, collect(zip(rep_strains, rep_profiles)))

        #/ Check if donor x recipient combination result in valid transconjugant(s), and if the recipient exist in current system
        for (recipient_index_Γ, (_,_)) in enumerate(zip(rep_strains, rep_profiles))
            #/ Check if the recipient exist in current system, and if donor x recipient combination result in valid transconjugant(s)
            j = findfirst(xy -> first(xy) == rep_strains[recipient_index_Γ] && last(xy) == rep_profiles[recipient_index_Γ], collect(zip(subpopulations, plasmid_profiles)))
            if any(propensities[:,donor_index_Γ,recipient_index_Γ] .> 0) && j !== nothing
                push!(recipients_of_donor, j)
                push!(recipient_indices_Γ, recipient_index_Γ)
            end    
        end      

        push!(recipients, recipients_of_donor)
        push!(recipients_indices_Γ, recipient_indices_Γ)
    end
    
    #/ Create and return struct
    _states = states(
        N, abundances, subpopulations,
        all_profiles, plasmid_profiles, donors, donors_indices_Γ, recipients, recipients_indices_Γ, rep_strains, rep_profiles
    )
    return _states 
end

"Wrapper function that gets all the per-capita rates of all the possible elementary reactions"
function get_rates(time, states, params)
    #/ Get background rates
    growth_rates, death_rates, carrying_capacities, segregation_rates = get_background_rates(time, 
		states, params.growth_rate, params.death_rate, params.perturbation_impact, params.carrying_capacity,
		params.segregation_error, params.plasmid_resistance, params.plasmid_cost
    )
    #/ Get interaction rates
    competition_rates, infection_rates = get_interaction_rates2(states, params)
    #/ Construct the struct and return it
    _rates = rates(
        growth_rates,
        death_rates,
        segregation_rates,
        competition_rates,
        infection_rates
    )
    return _rates 
end

#/ As the rates in the parameter file are defined per strain, we need to
#  broadcast the rates to their appropriate lengths here.
function get_background_rates(
    time,
    states,
    growth_rate,
    death_rate,
    perturbation,
    carrying_capacity,
    segregation_error,
    plasmid_resistance,
    plasmid_cost
    )
    #/ Allocate arrays of type Float64
    growth_rates = Array{Float64}(undef, 0)
    death_rates = Array{Float64}(undef, 0)
    carrying_capacities = Array{Float64}(undef, 0)
    segregation_rates = Array{Float64}(undef, 0)
	for (strain_id, plasmid_profile) in zip(states.subpopulations, states.plasmid_profiles)
        #/ Compute effective strain- and plasmid (infection state)-specific rates for each subpopulation
        #~ Growth rate: not including intra-specific competition 
        _plasmids = findall(plasmid_profile .> 0) # list the plasmids hosted by the focal subpopulation 
        
        η = (isempty(_plasmids)) ?
        growth_rate[strain_id] :
        growth_rate[strain_id] * reduce((x, y) -> x * (1.0 - y), plasmid_cost[_plasmids], init=1.0) # taking the sequential product of (1-cost) from each plasmid

        push!(growth_rates, η)
        
        #~ Death rate
        # Exclude time-dependent perturbation and plasmid resistance:
        # μ = death_rate[strain_id]

        # Include time-dependent perturbation and plasmid resistance:
        if perturbation[strain_id] == 0
            μ = death_rate[strain_id] 
        else
            μ = (isempty(_plasmids)) ? 
            death_rate[strain_id] * (1 + pulse_nd_single(time, perturbation[strain_id], 48)) :
            death_rate[strain_id] * (1 + pulse_nd_single(time, perturbation[strain_id], 48) * (1 - maximum(plasmid_resistance[_plasmids])))    
        end

        push!(death_rates, μ)
        #~ Carrying capacity
        C = carrying_capacity[strain_id]
        push!(carrying_capacities, C)
        #~ Segregation error
        ω = length(_plasmids) > 0 ? segregation_error[strain_id] : 0
        push!(segregation_rates, ω*η)
    end
    #/ Make the rate struct and return    
    return growth_rates, death_rates, carrying_capacities, segregation_rates
end

## Function that imposes time-dependent perturbation impact
#/ single pulse + no decrease, with perturbation parameter (treated as the fold of μ_i as strenth of pulse), and pulse time t_p (hour) defined
function pulse_nd_single(time, perturb, t_p)
    if time < t_p
       return 0.0
    else
       return perturb     
    end
end

## Function to update the death rates with time-dependent perturbatoin impact
function update_background_death(
    time,
    states,
    death_rate,
    perturbation,
    plasmid_resistance
    )
    #/ Allocate arrays of type Float64
    death_rates = Array{Float64}(undef, 0)
	for (strain_id, plasmid_profile) in zip(states.subpopulations, states.plasmid_profiles)
        #/ Compute effective strain- and plasmid (infection state)-specific rates for each subpopulation
        _plasmids = findall(plasmid_profile .> 0)
        #~ Death rate
        # Include time-dependent perturbation and plasmid resistance:
        if perturbation[strain_id] == 0
            μ = death_rate[strain_id]
        else    
            μ = (isempty(_plasmids)) ? 
            death_rate[strain_id] * (1 + pulse_nd_single(time, perturbation[strain_id], 48)) :
            death_rate[strain_id] * (1 + pulse_nd_single(time, perturbation[strain_id], 48) * (1 - maximum(plasmid_resistance[_plasmids])))
        end
            
        push!(death_rates, μ)
    end
    #/ Make the rate struct and return    
    return death_rates
end

## Function to update the per-capita infection rate (used for the new equation, as it will be used later to sample the donor in infection events)
function update_infection_rate(states, params, new_abundances)
    infection_rates = Array{Float64, 1}(undef, 0)
    for  (focal_substrain_idx,(focal_substrain, focal_plasmid_profile)) in enumerate(zip(states.subpopulations, states.plasmid_profiles))
        if any(focal_plasmid_profile .> 0)
            # note: If the focal subpouplation is a donor (i.e. has at least 1 plasmid), infection rate (γ_{i,p}) = strain-specific values, otherwise γ_{i,p} = 0
            # note: focal_sbustrain_idx is the subpopulation index in the system, where findfirst(x -> x == focal_substrain_idx, states.donors) is the subpopulation index in the donor list
            _infection_rate = params.infection_rate[focal_substrain] * sum(new_abundances[reduce(vcat, states.recipients[findfirst(x -> x == focal_substrain_idx, states.donors)])]) 
            #/ Push per-capita infection rate
            push!(infection_rates, _infection_rate) # per capita infection rate
        else
            #/ Otherwise, no infect, so push 0
            push!(infection_rates, 0.0)
        end
    end    
    return infection_rates
end    

"Compute the relevant interaction rates from the matrices/tensors"
function get_interaction_rates2(states, params)
    #/ Instantiate per capita competition rate and infection rate
    competition_rates = Array{Array{Float64, 1}, 1}(undef, 0)
    # infection_rates = Array{Array{Float64, 1}, 1}(undef, 0) # for the old infection equation, infection_rates is a vector of vectors
    infection_rates = Array{Float64, 1}(undef, 0) # for the new infection equation, infection_rates is a vector

    #/ Generate this K & matrix B based on matrix A, where B_ij = A_ij + 1 when i != j #(required when considering community-wise K instead of K_i)
    K = params.carrying_capacity[1]
    B = copy(params.A) .+1
    B[diagind(B)] .= 1
    

    
    #/ Generate an array that for each substrain contains the array of competition & infection rates
    # note: the density-dependent terms in the equations of per capita competition rate is taken care of later while calculating the total rates (see gillepie.jl; search for "Compute rates")
    for  (focal_substrain_idx,(focal_substrain, focal_plasmid_profile)) in enumerate(zip(states.subpopulations, states.plasmid_profiles)) 
        #/ Allocate arrays
        _competition_rates = Array{Float64, 1}(undef, 0)
    
        #/ Find those that carry at least one plasmid
        _plasmids = findall(focal_plasmid_profile .> 0) # return the indices (ids) of plasmid(s) the focal substrain hosts
        strain_id = focal_substrain[begin]        
        for other_substrain in states.subpopulations
            competitor_id = other_substrain[begin] #/ get the competitior's (or the non-focal substrain's) strain id
            #/ Compute the per-capita competition rate divided by carrying capacity: condider both intraspecific and interspecific competition and K (e.g. = K_i)
            η = (isempty(_plasmids)) ?
            params.growth_rate[strain_id] :
            params.growth_rate[strain_id] * reduce((x, y) -> x * (1.0 - y), params.plasmid_cost[_plasmids], init=1.0) # taking the sequential product of (1-cost) from each plasmid           

            _A = η * B[strain_id,competitor_id] / K
            
            push!(_competition_rates, _A)
        end
        push!(competition_rates, _competition_rates)

        #/ Calculate per-capita infection rate, i.e γ_{i,p} * sum (recipient abundance)
        if any(focal_plasmid_profile .> 0)
            # note: If the focal subpouplation is a donor (i.e. has at least 1 plasmid), infection rate (γ_{i,p}) = strain-specific values, otherwise γ_{i,p} = 0
            # note: focal_sbustrain_idx is the subpopulation index in the system, where findfirst(x -> x == focal_substrain_idx, states.donors) is the subpopulation index in the donor list
            _infection_rate = params.infection_rate[focal_substrain] * sum(states.abundances[reduce(vcat, states.recipients[findfirst(x -> x == focal_substrain_idx, states.donors)])]) 
            #/ Push per-capita infection rate
            push!(infection_rates, _infection_rate) # per capita infection rate
        else
            #/ Otherwise, no infect, so push 0
            push!(infection_rates, 0.0)
        end
    end

    #/ Return interaction rates
    return competition_rates, infection_rates
end

"Compute propensity tensor Γ"
# The propensity tensor P is used to show the propensities of recipient subpopulations, 
# being infected by the donor subpouplations, becoming transconjugant subpouplations
# The 3-dimention propensity tensor has the size = nm * nm * nm where n = number of all possilbe plasmid profiles and m = number of all possible bacterial strains.
# Note: this function used sparse array to save memory and increase iteration efficiency. 
function get_propensities3(params; power::Float64 = 2.0)
    #/ To compute the propensity tensor P, we need a list of all different strain ID & plasmid
    #  profiles that can exist in the system, even those that are not initially present.
    #!Note: The order of plasmid profiles is determined by Iterators.product.
    #!Note: Right now we can only define the propensity tensor here, but we would like it
    #       to be an input parameter to the system instead.
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
    end_time = now()
    elapsed_seconds = Dates.value(end_time - start_time)/1000.0
    @info "time elasped for this tensor generation (s):" elapsed_seconds 
    return Γ
end

end # module System1
#/ End module
