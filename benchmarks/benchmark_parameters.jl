const SEED    = 5357
const MODELS  = ["kendall", "mmek", "autoreg", "dimer-decay", "yeast"]

const ParamSet = Tuple{Float64,Int,Int,Int,Int}

parameters = Dict{String,ParamSet}()

##### Model 1 #####

# simulation parameters
t_final  = 4.0
n_saves  = 1000
n_trial  = 10

# benchmark parameters
n_sample = 1000
t_limit  = 30 * 60

parameters[MODELS[1]] = (t_final, n_saves, n_trial, n_sample, t_limit)

##### Model 2 #####

# simulation parameters
t_final  = 50.0
n_saves  = 1000
n_trial  = 10

# benchmark parameters
n_sample = 1000
t_limit  = 30 * 60

parameters[MODELS[2]] = (t_final, n_saves, n_trial, n_sample, t_limit)

##### Model 3 #####

# simulation parameters
t_final  = 500.0
n_saves  = 1000
n_trial  = 10

# benchmark parameters
n_sample = 1000
t_limit  = 30 * 60

parameters[MODELS[3]] = (t_final, n_saves, n_trial, n_sample, t_limit)

##### Model 4 #####

# simulation parameters
t_final  = 100.0
n_saves  = 1000
n_trial  = 10

# benchmark parameters
n_sample = 1000
t_limit  = 30 * 60

parameters[MODELS[4]] = (t_final, n_saves, n_trial, n_sample, t_limit)

##### Model 5 #####

# simulation parameters
t_final  = 100.0
n_saves  = 1000
n_trial  = 10

# benchmark parameters
n_sample = 1000
t_limit  = 30 * 60

parameters[MODELS[5]] = (t_final, n_saves, n_trial, n_sample, t_limit)
