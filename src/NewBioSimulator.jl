module NewBioSimulator

using Random, DataStructures
using StaticArrays

## imports

import Combinatorics: multiexponents
import MacroTools: postwalk
import StatsBase: sample
import Base: show, <=
import Base: getindex, setindex!, iterate, firstindex, lastindex, length, isempty, haskey, get
import Base: eachindex, push!, empty!, copy
import Base.Order: Ordering, ForwardOrdering, ReverseOrdering, Forward, Reverse

const global NBTYPES = [:vonneumann, :hexagon]

## load order

## interface

#### well-mixed
include(joinpath("interface", "well-mixed", "species.jl"))
include(joinpath("interface", "well-mixed", "reaction.jl"))
include(joinpath("interface", "well-mixed", "network.jl"))

export Species, Reaction, Network
export <=, number_species, number_reactions, species_list, reaction_list
export get_species, get_reaction

#### discrete-spatial
include(joinpath("interface", "spatial-discrete", "reaction.jl"))

export @def_reactions, @enumerate_with_sclass, @enumerate_with_nclass

## data structures
include(joinpath("data-structures", "dep_graph.jl"))

export DGView, DGLazy, DGVector, dependents, rxnrxn_depgraph, spcrxn_depgraph

include(joinpath("data-structures", "priority_queue.jl"))

export PQBinaryHeap, peektop

## state
include(joinpath("state", "site.jl"))
include(joinpath("state", "neighborhood.jl"))
include(joinpath("state", "vonneumann.jl"))
include(joinpath("state", "sample_class.jl"))
include(joinpath("state", "neighborhood_class.jl"))
include(joinpath("state", "hexagonal.jl"))
include(joinpath("state", "abstract_lattice.jl"))
include(joinpath("state", "SLattice.jl"))
include(joinpath("state", "NLattice.jl"))

export VonNeumann1D, VonNeumann2D, VonNeumann3D, Hexagonal2D
export SLattice, NLattice

## model
include(joinpath("model", "reaction_system.jl"))
include(joinpath("model", "ips.jl"))

## algorithms
include(joinpath("algorithms", "ssa_utils.jl"))
include(joinpath("algorithms", "abstract_algorithms.jl"))
include(joinpath("algorithms", "direct.jl"))
include(joinpath("algorithms", "enhanced_direct.jl"))
include(joinpath("algorithms", "sorting_direct.jl"))
include(joinpath("algorithms", "firstreaction.jl"))
include(joinpath("algorithms", "nextreaction.jl"))
include(joinpath("algorithms", "rejection.jl"))

export HasRates, HasSums

## simulators
include(joinpath("simulators", "abstract_simulators.jl"))
include(joinpath("simulators", "exact.jl"))
include(joinpath("simulators", "build_simulator.jl"))

export Direct, EnhancedDirect, SortingDirect
export FirstReaction, NextReaction, RejectionSSA

include("simulate.jl")

export simulate

end # module
