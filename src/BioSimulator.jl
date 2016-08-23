module BioSimulator

using Distributions

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

# Backend
include(joinpath("backend","reactionsystem.jl"))
include(joinpath("backend","pvec.jl"))
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

# Output
include(joinpath("output","partial_history.jl"))

# Simulate
include(joinpath("interface","simulate.jl"))

export Species, Reaction, Parameter, Network, simulate, petrinet, n_species, n_reactions, n_parameters, species_list, reaction_list, parameter_list, SSA, FRM, NRM, ODM, SAL

end # module
