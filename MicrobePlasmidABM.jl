#= MicrobePlasmid ABM =#
#~ Agent based model of the microbe-plasmid (host-pathogen) system, using community-wise carrying capacity
#~ Code used for HFSP project
#/ Start module
module MicrobePlasmidABM

#/ Local modules
#~ Parameter related modules
include("parameters/params.jl")
#~ Gillespie algorithm related modules
include("gillespie/system.jl")
include("gillespie/gillespie.jl")
include("gillespie/system1.jl")
include("gillespie/system2.jl")
include("gillespie/system3.jl")
include("gillespie/gillespie1.jl")
include("gillespie/gillespie2.jl")
include("gillespie/gillespie3.jl")
#~ Model related modules
#include("models/models.jl")
#~ Data I/O related modules
include("analysis/gatheroutput.jl")
include("analysis/gatheroutput2.jl")
#~ JSON file generation modules
include("parameters/json/json_generator.jl")
#~ Other modules
include("utils/utils.jl")

#/ Export some modules
#!Note: it might be more clear to omit these exports and always call
#       MicrobePlasmidABM.<ModuleName>
export Params
export Models

export Utils 
export DataHandler


# "include" the julia files for intellisense to work
macro ignore(args...) end
@ignore include("../scripts/examples.jl")

end # module MicrobePlasmidABM
#/ End module
