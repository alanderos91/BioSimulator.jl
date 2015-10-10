module BioSimulator

using Distributions
using DataFrames
using Gadfly
using LightGraphs

# Model Types
include(joinpath("model","species.jl"))
include(joinpath("model","reaction.jl"))
include(joinpath("model","network.jl"))
include(joinpath("model","simulation.jl"))

# Kinetics
include(joinpath("kinetics","mass_action.jl"))

# Output
include(joinpath("output","pstate.jl"))
include(joinpath("output","ptrajectory.jl"))
include(joinpath("output","simresult.jl"))
include(joinpath("output","simjob.jl"))
include(joinpath("output","util.jl"))
#include(joinpath("output","plot.jl"))

# Interface
include("interface.jl")

# Algorithms
include(joinpath("algorithm","ssa.jl"))
include(joinpath("algorithm","ossa.jl"))
include(joinpath("algorithm","frm.jl"))
include(joinpath("algorithm","nrm.jl"))
include(joinpath("algorithm","sal.jl"))

end # module
