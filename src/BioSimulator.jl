module BioSimulator

__precompile__()

using DataStructures
using RecipesBase
using SparseArrays
using LightGraphs: DiGraph, add_edge!
using Random

import TikzGraphs
import StatsFuns.RFunctions: poisrand
import DataFrames: DataFrame
import Statistics: mean, std

import Base: (<=),
             (>=),
             show,
             getindex,
             setindex!,
             isempty,
             length,
             size,
             eachindex,
             enumerate,
             fill!

# Backend
include(joinpath("backend","pvec.jl"))
include(joinpath("backend","reactionsystem.jl"))
include(joinpath("backend","dense.jl"))
include(joinpath("backend","sparse.jl"))

# Algorithms
include(joinpath("algorithm","algorithm.jl"))
include(joinpath("algorithm","ssa.jl"))
include(joinpath("algorithm","odm.jl"))
include(joinpath("algorithm","frm.jl"))
include(joinpath("algorithm","nrm.jl"))
include(joinpath("algorithm","otl.jl"))
include(joinpath("algorithm","sal.jl"))

# Interface
include(joinpath("interface","build_algorithm.jl"))
include(joinpath("interface","species.jl"))
include(joinpath("interface","reaction.jl"))
include(joinpath("interface","network.jl"))
include(joinpath("backend","util.jl"))

export Direct, FirstReaction, NextReaction, OptimizedDirect, TauLeaping, StepAnticipation
export Species, Reaction, Network, n_species, n_reactions, species_list, reaction_list

# Output
include(joinpath("output","petrinet.jl"))
include(joinpath("output","sample_path.jl"))
include(joinpath("output","regular_path.jl"))
include(joinpath("output","average_path.jl"))
include(joinpath("output","histogram.jl"))
include(joinpath("output","summary.jl"))

export visualize

include(joinpath("interface","simulate.jl"))

export simulate

end # module
