module BioSimulator

using Random, DataStructures, Statistics, PoissonRandom
using StaticArrays, SparseArrays, LinearAlgebra
using RecipesBase

## imports
import RecursiveArrayTools: DiffEqArray

import Combinatorics: multiexponents
import MacroTools: postwalk
import StatsBase: sample
import Base: show, <=
import Base: getindex, setindex!, iterate, firstindex, lastindex, length, isempty, haskey, get, in
import Base: eachindex, push!, empty!, copy
import Base.Order: Ordering, ForwardOrdering, ReverseOrdering, Forward, Reverse

## load order

## interface
include(joinpath("interface", "kineticlaw.jl"))

#### well-mixed
include(joinpath("interface", "well-mixed", "species.jl"))
include(joinpath("interface", "well-mixed", "reaction.jl"))
include(joinpath("interface", "well-mixed", "network.jl"))

export Species, Reaction, Network
export <=, number_species, number_reactions, species_list, reaction_list
export get_species, get_reaction

#### discrete-spatial
include(joinpath("interface", "spatial-discrete", "reaction.jl"))

# export @def_reactions, @enumerate_with_sclass, @enumerate_with_nclass
export @def_reactions, @enumerate_with_sclass

## data structures
include(joinpath("data-structures", "dep_graph.jl"))

export DGView, DGLazy, DGVector, dependents, rxnrxn_depgraph, spcrxn_depgraph

include(joinpath("data-structures", "priority_queue.jl"))

export PQBinaryHeap, peektop

## state
include(joinpath("state", "site.jl"))
include(joinpath("state", "vonneumann.jl"))
include(joinpath("state", "hexagonal.jl"))
include(joinpath("state", "Lattice.jl"))

export VonNeumann, Hexagonal
export Lattice

const global NBTYPES = [VonNeumann(), Hexagonal()]

## model
include(joinpath("model", "reaction_system.jl"))
include(joinpath("model", "SampleClassEnumeration.jl"))
include(joinpath("model", "ips.jl"))

## algorithms
include(joinpath("algorithms", "ssa_utils.jl"))
include(joinpath("algorithms", "abstract_algorithm.jl"))
include(joinpath("algorithms", "direct.jl"))
include(joinpath("algorithms", "enhanced_direct.jl"))
include(joinpath("algorithms", "sorting_direct.jl"))
include(joinpath("algorithms", "firstreaction.jl"))
include(joinpath("algorithms", "nextreaction.jl"))
include(joinpath("algorithms", "rejection.jl"))

include(joinpath("algorithms", "tauleap_utils.jl"))
include(joinpath("algorithms", "poisson_tauleap.jl"))
include(joinpath("algorithms", "stepanticipation.jl"))

export HasRates, HasSums

## simulators
include(joinpath("simulators", "abstract_simulator.jl"))
include(joinpath("simulators", "exact.jl"))
include(joinpath("simulators", "tauleap.jl"))
include(joinpath("simulators", "tauleap_hybrid.jl"))
include(joinpath("simulators", "build_simulator.jl"))

export Direct, EnhancedDirect, SortingDirect
export FirstReaction, NextReaction, RejectionSSA
export TauLeapingDG2001, TauLeapingDGLP2003, StepAnticipation
export HybridSAL

export cumulative_intensity, jump_rates, next_jump_index, next_jump_time

## simulation + output
include("simulate.jl")
include(joinpath("output", "SamplePath.jl"))
include(joinpath("output", "Configuration.jl"))
include(joinpath("output", "tabletraits.jl"))

export simulate, parse_model, SamplePath, Ensemble, Configuration
export save_state, save_rates

end # module
