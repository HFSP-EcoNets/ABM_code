#= Utility functions =#

#/ Start module
module Utils

#/ Packages
using CSV
using DataFrames
using ProgressBars
#using WeakRefStrings

#/ Functions

# Functions that replaces a string value ("NA") in column "job_id" by another string value in a csv file 
# Option 1: using DataFrames
function create_append_job_id(job_id_path::String, job_key::String, job_id::String)
    # Create a new .csv file with the job_id_path, and append the job_key and job_id
        CSV.write(job_id_path, DataFrame(key = job_key, job_id = job_id), append=true)   
end

end # module Utils
#/ End module


# Function that demonstrate progress bar with a given threshold of iteration (random, float intervals)
function progressbar(t_final)
    
total_iterations = t_final

pb = ProgressBar(total=total_iterations)

iteration = 0

while iteration < total_iterations

    del_t = rand()  
    
    sleep(0.1)
    if (iteration+del_t > ceil(iteration) )
        step = Int(floor(iteration+del_t))-Int(floor(iteration))
        update(pb,step)
        
    end
    iteration += del_t
    
end

end