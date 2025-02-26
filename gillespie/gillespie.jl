#= Implementation of the Gillespie algorithm for agent-based model, using community-wise carrying capacity K and competition matrix B =#
#/ Start module
module Gillespie

#/ Packages
using Random
# using Distributions
using StatsBase
using Dates
using DataFrames
using ProgressBars
using LinearAlgebra
# using SQLite: DB, Stmt, bind!, load!
# using SQLite.DBInterface: execute
#/ Local packages
using MicrobePlasmidABM

### FUNCTIONS
function run(
    params::PType;
    filename::String,
    output_folder::String,
    job_key::String,
    JOB_ID::String
    ) where {PType <: MicrobePlasmidABM.Params.params}
    start_time = now()     #/ Record start time

    #/ Set random seed
    seed = (params.rng_seed === nothing) ? 42 : params.rng_seed
    # Random.seed!(seed) # Approach of Han
    rng = MersenneTwister(seed)
    #/ Set some parameters
	  t = 0.0
    twrite = params.t_output

    pb = (Sys.isapple()) ? ProgressBar(total=Int(params.t_final)) : ProgressBar(total=Int(params.t_final), output_stream=stdout) # the later is for progress bar in HPC; Sys.islinux()

    #/ Initialize system
    #~ Get the states and rates structs, and the propensity StaticArray (tensor)
    @info "Initializing system..."
    states, rates, propensities = MicrobePlasmidABM.System.initialize_system(t, params)    
    
    #/ Make array that will always hold the change in abundances, as it will be needed
    #  to know how the rates should be updated.
    new_abundances = copy(states.abundances)
    
    #/ Compute rates
    @info "Computing rates"
    Rgrowth = sum(rates.growth_rates .* states.abundances)
    Rdeath = sum(rates.death_rates .* states.abundances)
    Rsegregation = sum(rates.segregation_rates .* states.abundances)
    Rcompetition = sum(
        [sum(rates .* states.abundances) for rates in rates.competition_rates] .*
            states.abundances
    )
    
    Rinfection = sum(rates.infection_rates .* states.abundances)    


    #/ Compute total rate
    Rtotal = Rgrowth + Rdeath + Rcompetition + Rsegregation + Rinfection

    #/ Define function array
    events = [
        do_growth!,
        do_death!,
        do_segregation!,
        do_competition!,
        do_infection!
    ]
    n_events = zeros(Int, length(events))

    #/ Create an empty DataFrame to store
    # (1) the system state (abundance)
    # (2) the event counts at desired time step
    state_df = DataFrame(
        t = repeat([t], states.N), #/time, repeated with number of subpopulations
        subpop_id = 1:states.N, #/subpouplation ids
        strain = states.subpopulations, #/strains of subpopulaitons
        p_profile = [string(v) for v in states.plasmid_profiles], #/plasmid profiles of subpoupations (in string type) 
        abundance = states.abundances #/abundance of subpopulations
    )
    event_df = DataFrame(
        t = t,
        growth = n_events[1],
        death = n_events[2],
        segregation = n_events[3],
        competition = n_events[4],
        infection = n_events[5]
    )

    #/ Initialize database (don't print the returned variable)
    @info "Initializing database"
    MicrobePlasmidABM.DataHandler.initialize_database(output_folder, job_key);
    
    #/ Execute Gillespie algorithm until time t_final
    @info "Begin simulation loop"

    while t < params.t_final
        #/ Determine time step from exponential distribution
        Δt = randexp(rng) / Rtotal # approach of Armun (note: randexp() only allows rate = 1 for the exponential distribution)

        #/ Sample state change using StatsBase
        #~ Sample a specific event weighted by their rate proportion
        _event_idx = sample(rng,
            1:length(events),
            Weights([Rgrowth,Rdeath,Rsegregation,Rcompetition,Rinfection])
        )
        #/ Execute stage change
        n_events[_event_idx] += 1
        #!Note: Infection is the only function that needs the propensities, so check here
        #       if infection is chosen and change the function call accordingly.
        #!Note: The check assumes that infection is the 5th event.
        affected_subpopulations = (_event_idx == 5) ?
            events[_event_idx](rng, states, rates, new_abundances, propensities, params, t) : #/ execute do_infection!(states, rates, new_abundances, propensities, params, t); return recipient's and transconjugant's id
            events[_event_idx](rng, states, rates, new_abundances) #/ execute do_growth/death/segregation/competition!(states, rates, new_abundances); return focal subpopulation's id
        #/ Update rates
        Rgrowth, Rdeath, Rsegregation, Rcompetition, Rinfection = update_rates!(
            t,
            states, 
            params,
            rates, 
            affected_subpopulations,
            new_abundances,
            Rgrowth,
            Rdeath,
            Rsegregation,
            Rcompetition,
            Rinfection
        )
        Rtotal = Rgrowth + Rdeath + Rsegregation + Rcompetition + Rinfection

        #/ Udate progress bar
        if (t+Δt > ceil(t) )
            step = Int(floor(t+Δt))-Int(floor(t))
            update(pb,step)
            
        end

        #/ Increment time
        t = t + Δt

        #/ Store desired data
        if t > twrite
            

            # Store data into dataframe
            MicrobePlasmidABM.DataHandler.write_current_state(t, states, state_df)
            MicrobePlasmidABM.DataHandler.write_current_event(t, n_events, event_df)

            # Store dataframe to SQLite database every n = t_final steps
            if t % params.t_final < params.t_output
                MicrobePlasmidABM.DataHandler.save_df_to_sqlite(state_df, event_df, output_folder, job_key)
            end

            twrite += params.t_output
        end

        #/ Break & store the dataframe to SQLite database if abundances are 0
        if sum(states.abundances) == 0
            MicrobePlasmidABM.DataHandler.save_df_to_sqlite(state_df, event_df, output_folder, job_key)
            @info "Full extinction, stopping early."
            break
        end

        flush(stdout) #a function call that flushes the output buffer associated with the standard output stream 'stdout'.
    end

    #/ Write data (seed, filename, start_time, end_time, elapsed_seconds)
    end_time = now()
    elapsed_seconds = Dates.value(end_time - start_time)/1000.0
    MicrobePlasmidABM.DataHandler.save_info_to_sqlite(
        seed, filename, start_time, end_time, elapsed_seconds, output_folder, job_key, JOB_ID
    )

end


"Update the rates"
function update_rates!(
    time, states, params, rates, 
    affected_subpopulations,
    new_abundances,
    Rgrowth,
    Rdeath,
    Rsegregation,
    Rcompetition,
    Rinfection
    )
    #/ Compute the effective change in abundances
    Δabundance = new_abundances .- states.abundances

    #/ Update background (per capita) death rate
    rates.death_rates = MicrobePlasmidABM.System.update_background_death(
        time, 
        states,
        params.death_rate,
        params.perturbation_impact,
        params.plasmid_resistance)      

    #/ Loop over all the affected subpopulations and update the total rates accordingly
    #!Note: most elementary reactions result in the change of a single subpopulation,
    #       but infection give rise to changes in multiple subpopulations.
    #       This is the reason for the loop.
    for k in affected_subpopulations
        #/ Update background rates
        Rgrowth = Rgrowth + Δabundance[k] * rates.growth_rates[k]
        Rdeath = sum(rates.death_rates .* states.abundances)  + Δabundance[k] * rates.death_rates[k]
        Rsegregation = Rsegregation + Δabundance[k] * rates.segregation_rates[k]
        
        #/ Update interaction rates
        #~ Competition
        ΔRcompetition = sum(
            [new_abundances[j] .* rates.competition_rates[j][k] for j in 1:states.N]
        )
        ΔRcompetition += sum(states.abundances .* rates.competition_rates[k])
        Rcompetition = Rcompetition + Δabundance[k] * ΔRcompetition

        #~ Infection (considering recipient subpouplations in per-capita rate)
        # (update per-capita infection rates of affected subpopulations, then recalculate the total infection rate)
        for (i, (donor)) in enumerate(states.donors)
            if k in states.recipients[i] # change in donor's recipient's abundance (k as the focal subpopulation's recipient in the infection equation)
                Δinfection_rate = params.infection_rate[states.subpopulations[donor]] * Δabundance[k]
                rates.infection_rates[donor] += Δinfection_rate 
            end      
        end

        
    end
      
      Rinfection = sum(rates.infection_rates .* new_abundances) # this line is used for both Attempt 2 & 3

    #/ Update the abundance
    states.abundances .+= Δabundance
    #/ Check for rounding errors
    return Rgrowth, Rdeath, Rsegregation, Rcompetition, Rinfection
end

"Growth event function"
function do_growth!(rng, states, rates, new_abundances)       
    #/ Sample subpopulation k to increase
    k = sample(rng, 1:states.N, Weights(rates.growth_rates .* states.abundances))
    new_abundances[k] = new_abundances[k] + 1
    #/ Just return affected subpopulations
    return [k]
end

"Death event function"
function do_death!(rng, states, rates, new_abundances)
    #/ Sample population to decrease
    k = sample(rng, 1:states.N, Weights(rates.death_rates .* states.abundances))
    new_abundances[k] = max(0, new_abundances[k] - 1) 
    return [k]
end

"Competition event function"
function do_competition!(rng, states, rates, new_abundances)
    #/ Sample the subpopulation that will suffer from competition 
    _Rcompetition = [sum(states.abundances .* rates.competition_rates[i]) for i in 1:states.N] # get the latest per capita competition rate for each subpoupation in the system
    k = sample(rng, 1:states.N, Weights(_Rcompetition .* states.abundances))
    #/ Update abundances
    new_abundances[k] = max(0, new_abundances[k] - 1)
    return [k]
end

"Segregation event"
function do_segregation!(rng, states, rates, new_abundances)
    #/ Sample parent strain that will lose all plasmids
    parent = sample(rng, 1:states.N, Weights(rates.segregation_rates .* states.abundances))
    #/ Find the subpopulation of the strain that has no plasmids, if it exists
    _zero_profile = zeros(Int, length(states.plasmid_profiles[begin]))
    
    child = findfirst(
        states.subpopulations .== states.subpopulations[parent] .&&
            states.plasmid_profiles .== [_zero_profile]
    )
    if !isnothing(child)
        new_abundances[child] = new_abundances[child] + 1
        return [child]
    else
        #/ If it does not exist, create it and all other relevant elements
        child = states.N + 1
        _strain_id = states.subpopulations[parent][begin]
        _plasmid_profile = zeros(length(states.plasmid_profiles[parent]))

        #/ Update state arrays
        states.N = states.N + 1
        push!(states.abundances, 0)      # !Note: old abundances, so push 0
        push!(new_abundances, 1)         # New abundances
        push!(states.subpopulations, _strain_id)
        push!(states.plasmid_profiles, _plasmid_profile)
        #/ New subpopulation will always be a recipient as it has just lost all its plasmids
        for (i, _donor) in enumerate(states.donors)
            push!(states.recipients[i], states.N)
        end

        #/ Update rate arrays
        push!(rates.growth_rates, rates.growth_rates[parent])
        push!(rates.death_rates, rates.death_rates[parent])
        push!(rates.segregation_rates, 0.0)
        for i in 1:(states.N-1)
            push!(rates.competition_rates[i], rates.competition_rates[parent][i])
        end
        push!(rates.competition_rates, rates.competition_rates[parent])
    end
    return [child]
end

"Infection event"
function sampled_donor(rng, states, rates) # function that will return a sampled donor (index in the donor list) and its recipients (recipients' subpopulation id)
    while true
        donor_index = sample(rng, 1:length(states.donors), Weights(rates.infection_rates[states.donors] .* states.abundances[states.donors])) # get donors indice in the donors; for the new infection equation
        its_recipients = states.recipients[donor_index] #/ this return the element(s) (i.e. the recipients' subpopulation id)
        if any(x -> x > 0, states.abundances[its_recipients]) #/ if any of the sampled donor's recipients has abundance > 0; if no, resample the donor
            return donor_index, its_recipients
        end    
    end           
end    
function do_infection!(rng, states, rates, new_abundances, propensities, params, time)
    #/ Sample a donor and one of its recipients
    donor_id, recipients = sampled_donor(rng, states, rates)

    #~ Sample recipients based on their abundance /and/ infection rate
    if all(states.abundances[recipients] .== 0) # if all the recipients have become extinct
        # error("recipients are all extinct!") # report error 
    else       
        recipient = sample(rng,
            recipients,
            Weights(states.abundances[recipients])
        ) #/ this return the element(s) (i.e. the recipient's subpopulation id)
    
        # get recipient's indice in the recipients
        recipient_id = findfirst(x -> x == recipient, recipients)

        #~ Sample transconjugant based on the propensity tensor
        #!Note: We need to get the subpopulation index in propensity tensor,
        #       as the propensity tensor is nm * nm * nm where n = number of all possilbe plasmid profiles and m = number of all possible bacterial strains
        i = states.donors_indices_Γ[donor_id]
        j = states.recipients_indices_Γ[donor_id][recipient_id]
        _transconjugant = sample(rng, 1:params.n_bstrains * length(states.all_profiles), Weights(propensities[:, i, j])) # transconjugant's subpouplation indice in Γ    
        _strain = states.rep_strains[_transconjugant]
        _plasmid_profile = states.rep_profiles[_transconjugant]

        transconjugant = findfirst(xy -> first(xy) == _strain && last(xy) == _plasmid_profile, collect(zip(states.subpopulations, states.plasmid_profiles)))
        # @info transconjugant recipient donor # used to check if idices are reasonable 

        #/ Update abundances
        #~ Remove one of the recipient population
        new_abundances[recipient] = max(0, new_abundances[recipient] - 1)
        #/ Update state
        if !isnothing(transconjugant)
            #/ If it exists, update its abundance
            new_abundances[transconjugant] = new_abundances[transconjugant] + 1
            return [recipient, transconjugant]
        else
        #/ If it does not exist, create it and all other relevant elements
            transconjugant = states.N + 1
            _strain_id = states.subpopulations[recipient][begin]
            #/ Update state arrays
            states.N = states.N + 1
            push!(states.abundances, 0)      # !Note: old abundances, so push 0
            push!(new_abundances, 1)         # New abundances        
            push!(states.subpopulations, _strain_id)
            push!(states.plasmid_profiles, _plasmid_profile)


            # get the new subpouplation's indice in the propensity tensor (which is the same as a donor & as a recipient)
            _index_Γ = findfirst(xy -> first(xy) == _strain_id && last(xy) == _plasmid_profile, collect(zip(states.rep_strains, states.rep_profiles)))

        
            #/ New subpouplation could be a recipient of the pre-existing donors, so update thier recipient lists
            for i in states.donors
                donor_strain = states.subpopulations[i]
                donor_plasmid_profile = states.plasmid_profiles[i]
                donor_index_Γ = findfirst(xy -> first(xy) == donor_strain && last(xy) == donor_plasmid_profile, collect(zip(states.rep_strains, states.rep_profiles)))
                recipient_index_Γ = findfirst(xy -> first(xy) == _strain_id && last(xy) == _plasmid_profile, collect(zip(states.rep_strains, states.rep_profiles)))
                if any(propensities[:,donor_index_Γ,recipient_index_Γ] .> 0)
                    j = findfirst(x -> x == i, states.donors) # get the indice of the donor in the donor list
                    push!(states.recipients[j], transconjugant)
                    push!(states.recipients_indices_Γ[j], _index_Γ)
                end    
            end     

            #/ New subpopulation will always be a donor, as it has just gained a plasmid
            push!(states.donors, transconjugant)
            push!(states.donors_indices_Γ, _index_Γ)

            #~ So we gather its recipients here
            recipients_of_donor = Array{Int, 1}(undef, 0)
            _recipient_indices_Γ = Array{Int, 1}(undef, 0)

            for (recipient_index_Γ, (_,_)) in enumerate(zip(states.rep_strains, states.rep_profiles))
             #/ Check if the recipient exist in current system, and if donor x recipient combination result in valid transconjugant(s)
             j = findfirst(xy -> first(xy) == states.rep_strains[recipient_index_Γ] && last(xy) == states.rep_profiles[recipient_index_Γ], collect(zip(states.subpopulations, states.plasmid_profiles)))
             if any(propensities[:,_index_Γ,recipient_index_Γ] .> 0) && j !== nothing
                push!(recipients_of_donor, j)
                push!(_recipient_indices_Γ, recipient_index_Γ)
             end    
            end 

            push!(states.recipients, recipients_of_donor)
            push!(states.recipients_indices_Γ, _recipient_indices_Γ)

            #/ Update rate arrays
            _plasmids = findall(_plasmid_profile .> 0) # list the plasmids hosted by the new subpopulation

            # growth rate: not including intra-specific competition 
            _η = params.growth_rate[_strain_id] * reduce((x, y) -> x * (1.0 - y), params.plasmid_cost[_plasmids], init=1.0)
            push!(rates.growth_rates, _η)
        
            # death rate
            _μ = params.death_rate[_strain_id] * (1 + pulse_expd_periodic2(100.0, time, params.perturbation_impact[_strain_id], 100, 9.0) * (1 - maximum(params.plasmid_resistance[_plasmids])))
            push!(rates.death_rates, _μ)
	    # segregation rate	
            _ω = params.segregation_error[_strain_id]
            push!(rates.segregation_rates, _ω * _η)


            #/ competition rate: (1) extend existing vectors of competition rates with the old subpopulations as the focal subpopulations
            for i in 1:(states.N-1)
                #/ feed in the new subpopulation's competiton impact on old (focal) subpopulations (η * A)
                push!(rates.competition_rates[i], rates.competition_rates[i][recipient])
            end

            #/ Generate this K & matrix B based on matrix A, where B_ij = A_ij + 1 when i != j #(required when considering community-wise K instead of K_i)
            K = params.carrying_capacity[1]
            B = copy(params.A) .+1
            B[diagind(B)] .= 1
            
            #/ competition rate: (2) create new vectors of competition rates with the new subpopulation as the focal subpopulation
            _competition_rates = Array{Float64, 1}(undef, 0)
            
            for i in 1:(states.N)
                # consider both intraspecific and interpecific competition rate, and community-wise K 
                push!(_competition_rates, _η * B[_strain_id, states.subpopulations[i]]/ K)
                    
            end

            #/ infection rate: calculate the per-capita infection rate of the new subpopulation as the focal subpopulation
            _infection_rate = params.infection_rate[_strain_id] * sum(states.abundances[recipients_of_donor])
            
            #/ add the new vectors of competition rates and infectoin rates to the existing vectors of vectors
            push!(rates.competition_rates, _competition_rates)
            push!(rates.infection_rates, _infection_rate) 
        
        end
    end    
    return [recipient, transconjugant]
end

function pulse_expd_periodic2(lag::Any, time::Any, perturb::Any, t_p::Any, k::Any)
    if time < lag
        return 0.0
    else        
        return k * exp(-(1/perturb) * (time % t_p))
    end     
end

end # module Gillespie
#/ End module
