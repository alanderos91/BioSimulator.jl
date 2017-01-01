__precompile__(true)

module BioSimulator

using Distributions
using Reexport

@reexport using Plots

using GraphViz

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

const ALGORITHMS = [:ODM, :SSA, :FRM, :NRM, :SAL]

# Interface
include(joinpath("interface","species.jl"))
include(joinpath("interface","reaction.jl"))
include(joinpath("interface","network.jl"))

export Species, Reaction, Network,
             n_species, n_reactions, species_list, reaction_list

# Backend
include(joinpath("backend","pvec.jl"))
include(joinpath("backend","reactionsystem.jl"))
include(joinpath("backend","util.jl"))
include(joinpath("backend","dense.jl"))
include(joinpath("backend","sparse.jl"))

# Algorithms
include(joinpath("algorithm","algorithm.jl"))
include(joinpath("algorithm","ssa.jl"))
include(joinpath("algorithm","odm.jl"))
include(joinpath("algorithm","frm.jl"))
include(joinpath("algorithm","nrm.jl"))
include(joinpath("algorithm","sal.jl"))

export SSA, FRM, NRM, ODM, SAL

# Output
include(joinpath("output","partial_history.jl"))
include(joinpath("output","plot.jl"))
include(joinpath("output","petrinet.jl"))

export petrinet, get_data, get_metadata, save_data

# Simulate
include(joinpath("interface","simulate.jl"))

export simulate

end # module
