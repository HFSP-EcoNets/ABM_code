#= Handle data from simulations to store into desired formats =#

#/ Start module
module DataHandler

#/ Packages
using Dates
using DataFrames
using JSON
using Parameters
using SQLite: DB, Stmt, bind!, load!
using SQLite.DBInterface: execute

### FUNCTIONS
"Write current state into a dataframe (without writing into db)"
function write_current_state(t, states, state_df)
    # extract new elements
    times = repeat([t], states.N)
    subpop_ids = 1:states.N
    strains = states.subpopulations
    profiles = [string(v) for v in states.plasmid_profiles]
    abundances = states.abundances

    # add new elements to each column
    state_df.t = append!(state_df.t, times)
    state_df.subpop_id = append!(state_df.subpop_id, subpop_ids)
    state_df.strain = append!(state_df.strain, strains)
    state_df.p_profile = append!(state_df.p_profile, profiles)
    state_df.abundance = append!(state_df.abundance, abundances)
end

"Write current event counts into a dataframe (without writing into db)"
function write_current_event(t, n_events, event_df)
    # add new elements to each column
    event_df.t = append!(event_df.t, t)
    event_df.growth = append!(event_df.growth, n_events[1])
    event_df.death = append!(event_df.death, n_events[2])
    event_df.segregation = append!(event_df.segregation, n_events[3])
    event_df.competition = append!(event_df.competition, n_events[4])
    event_df.infection = append!(event_df.infection, n_events[5])
end 

"Write dataframes into SQLite database"
function save_df_to_sqlite(state_df, event_df, output_folder::String  = "output", job_key::String  = "000")
    # link to database
    db = DB(joinpath(output_folder,"output$job_key.sqlite"))

    # insert dataframes into the database's tables
    #load!(db, "bsubabundance", state_df)
    #load!(db, "events", event_df)
    tablename = state_df |> load!(db, "bsubabundance")
    tablename = event_df |> load!(db, "events")
end

"Write simulation information into SQLite database"
function save_info_to_sqlite(seed::Int, json, start_time, end_time, elapsed_seconds, output_folder::String  = "output", job_key::String  = "000", JOB_ID::String = "123456")

    # link to database
    db = DB(joinpath(output_folder,"output$job_key.sqlite"))

    # write simulation information into table "meta"
    execute(db, "BEGIN TRANSACTION")

    execute(db,
            "INSERT INTO meta VALUES (?,?)",
            ["seed", seed]
            )

    execute(db,
            "INSERT INTO meta VALUES (?,?)",
            ["key", job_key]
            )        
    
    execute(db,
            "INSERT INTO meta VALUES (?,?)",
            ["job", JOB_ID]
            ) 

    execute(db,
            "INSERT INTO meta VALUES (?,?)",
            ["json", json]
            )
            
    execute(db,
            "INSERT INTO meta VALUES (?,?)",
            ["start_time", Dates.format(start_time, "yyyy-mm-ddTHH:MM:SS")]
            )

    execute(db,
            "INSERT INTO meta VALUES (?,?)",
            ["end_time", Dates.format(end_time, "yyyy-mm-ddTHH:MM:SS")]
            )

    execute(db,
            "INSERT INTO meta VALUES (?,?)",
            ["elapsed_seconds", elapsed_seconds]
            )
    execute(db, "COMMIT")


end    

"Initialize SQLite database"
function initialize_database(
    output_folder::String = "output",
    job_key::String = "000";
    outfilename::String = joinpath(output_folder,"output$job_key.sqlite")
        )
    if isfile(outfilename)
        rm(outfilename)
        # error("output.sqlite already exists; delete first")
    end
    db = DB(outfilename)

    # create table "meta" in database: e.g. seed, job key, job id, json file used, start time, end time, elapsed seconds
    execute(db, "CREATE TABLE meta (key, value);")

    # create table "bsubabundance" for exsisting bacterial substrain: time, substrain ID, strain_id, plasmid_profile, abundance
    execute(db, """
        CREATE TABLE bsubabundance (
            t REAL,
            subpop_id INTEGER,
            strain INTEGER,
            p_profile STRING,
            abundance INTEGER
        )
    """)

    # create table "events" for accumulated number of events: time, grwoth, death, segregation, competition
    #(To be added: infection)
    execute(db, """
        CREATE TABLE events (
            t REAL,
            growth INTEGER,
            death INTEGER,
            segregation INTEGER,
            competition INTEGER,
            infection INTEGER
        )
    """)

    return db	  
end

end # module DataHandler
#/ End module
