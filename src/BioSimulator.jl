module BioSimulator

using Distributions
using DataFrames
using Gadfly
using LightGraphs

abstract OutputType
immutable Explicit  <: OutputType end
immutable Uniform   <: OutputType end
immutable Mean      <: OutputType end
immutable Histogram <: OutputType end

abstract Algorithm

# Interface
include(joinpath("interface","species.jl"))
include(joinpath("interface","reaction.jl"))
include(joinpath("interface","parameter.jl"))
include(joinpath("interface","network.jl"))
include(joinpath("interface","simulate.jl"))

# Kinetics
include(joinpath("kinetics","mass_action.jl"))

# Algorithms
include(joinpath("algorithm","ssa.jl"))
include(joinpath("algorithm","odm.jl"))
include(joinpath("algorithm","frm.jl"))
include(joinpath("algorithm","nrm.jl"))
include(joinpath("algorithm","sal.jl"))

# Output
include(joinpath("output", "update.jl"))

export Network, Simulation, Explicit, Uniform, Mean, Histogram,
  simulate, Species, Reaction, parameter
end # module
