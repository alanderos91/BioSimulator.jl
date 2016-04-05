module BioSimulator

using Distributions
using DataFrames
using Gadfly
using LightGraphs
using GraphViz
using Compose

import Base: getindex,
             setindex!,
             isempty,
             length,
             size,
             start,
             done,
             next,
             eachindex,
             enumerate,
             fill!
import Compose: MeasureOrNumber

abstract OutputType
immutable Explicit <: OutputType end
immutable Uniform  <: OutputType end

const ALGORITHMS = [:SSA, :SAL]

# Interface
include(joinpath("interface","species.jl"))
include(joinpath("interface","reaction.jl"))
include(joinpath("interface","parameter.jl"))
include(joinpath("interface","network.jl"))

# Backend
include(joinpath("backend","mass_action.jl"))
include(joinpath("backend","rxnvec.jl"))
include(joinpath("backend","model.jl"))

# Algorithms
include(joinpath("algorithm","algorithm.jl"))
include(joinpath("algorithm","ssa.jl"))
#include(joinpath("algorithm","odm.jl"))
#include(joinpath("algorithm","frm.jl"))
#include(joinpath("algorithm","nrm.jl"))
include(joinpath("algorithm","sal.jl"))

# Output
include(joinpath("output", "observer.jl"))
include(joinpath("output", "update.jl"))
include(joinpath("output", "util.jl"))
include(joinpath("output", "plot.jl"))
include(joinpath("interface", "petrinet.jl"))

include(joinpath("interface","simulate.jl"))

export Network, Simulation, Explicit, Uniform,
  simulate, Species, Reaction, parameter,
  petrinet, plot, species
end # module
