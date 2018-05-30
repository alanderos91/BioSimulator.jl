module BioSimulator

using DataStructures
using Reexport
using Interact
using Reactive

@reexport using Plots
@reexport using DataFrames

import TikzGraphs
using LightGraphs: DiGraph, add_edge!

import StatsFuns.RFunctions: poisrand

import Base: (<=),
             (>=),
             show,
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
include(joinpath("interface","gui.jl"))
include(joinpath("backend","util.jl"))

export Direct, FirstReaction, NextReaction, OptimizedDirect, TauLeaping, StepAnticipation
export Species, Reaction, Network,
             n_species, n_reactions, species_list, reaction_list

# export generate_gui, plot_interface

# Output
# include(joinpath("output","partial_history.jl"))
include(joinpath("output","petrinet.jl"))
# include(joinpath("output","stats.jl"))
# include(joinpath("output","plot.jl"))
include(joinpath("output","sample_path.jl"))
include(joinpath("output","regular_path.jl"))

# export get_data, get_dataframe, visualize
# export Trajectory, MeanTrajectory, Histogram, PhaseTrajectory

include(joinpath("interface","simulate.jl"))

export simulate

end # module
