module NewBioSimulator

using Random, DataStructures

## imports

import Base: show, <=
import Base: getindex, setindex!, iterate, firstindex, lastindex, length, isempty, haskey, get
import Base.Order: Ordering, ForwardOrdering, ReverseOrdering, Forward, Reverse

## load order

## interface
include(joinpath("interface", "well-mixed", "species.jl"))
include(joinpath("interface", "well-mixed", "reaction.jl"))
include(joinpath("interface", "well-mixed", "network.jl"))

export Species, Reaction, Network
export <=, number_species, number_reactions, species_list, reaction_list
export get_species, get_reaction

## data structures
include(joinpath("data-structures", "dep_graph.jl"))

export DGView, DGLazy, DGVector, dependents

include(joinpath("data-structures", "priority_queue.jl"))

export PQBinaryHeap, peektop

## model
include(joinpath("model", "reaction_system.jl"))

## algorithms
include(joinpath("algorithms", "ssa_utils.jl"))
include(joinpath("algorithms", "abstract_algorithms.jl"))
include(joinpath("algorithms", "direct.jl"))
include(joinpath("algorithms", "firstreaction.jl"))

export HasRates, HasSums

## simulators
include(joinpath("simulators", "abstract_simulators.jl"))
include(joinpath("simulators", "exact.jl"))
include(joinpath("simulators", "build_simulator.jl"))

export Direct, FirstReaction

include("simulate.jl")

export simulate

end # module
