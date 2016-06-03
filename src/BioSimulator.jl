module BioSimulator

using Gadfly
using Distributions
using LightGraphs
using DataFrames
using GraphViz

import Compose: cm, inch, mm, pt, px

import Base: (<=),
             (>=),
             getindex,
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

abstract OutputType
immutable Explicit <: OutputType end
immutable Uniform  <: OutputType end

const ALGORITHMS = [:ODM, :SSA, :FRM, :NRM, :SAL]

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
include(joinpath("algorithm","odm.jl"))
include(joinpath("algorithm","frm.jl"))
include(joinpath("algorithm","nrm.jl"))
include(joinpath("algorithm","sal.jl"))

# Output
include(joinpath("output", "observer.jl"))
include(joinpath("output", "update.jl"))
include(joinpath("output", "util.jl"))
include(joinpath("output", "plot.jl"))
include(joinpath("interface", "petrinet.jl"))

include(joinpath("interface","simulate.jl"))

export Species, Reaction, Parameter, Network, simulate, petrinet, n_species, n_reactions, n_parameters, species_list, reaction_list, parameter_list, value, id, get_speciesdf, get_metadata, plot_results

end # module
