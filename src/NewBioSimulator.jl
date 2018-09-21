module NewBioSimulator

using Random, DataStructures

## imports

import Base: show, <=

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

end # module
