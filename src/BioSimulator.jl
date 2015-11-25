module BioSimulator

using Distributions
using DataFrames
using Gadfly
using LightGraphs

abstract OutputType
immutable Explicit  <: OutputType end
immutable Uniform   <: OutputType end
immutable Mean      <: OutputType end
immutable Histogram <: OutputType end

abstract Algorithm

_init(alg::Algorithm, rxns, spcs, params) = return;
_reset(alg::Algorithm, rxns, spcs, params) = return;
_step(alg::Algorithm, rxns, spcs, params) = return;

# Interface
include(joinpath("interface","species.jl"))
include(joinpath("interface","reaction.jl"))
include(joinpath("interface","network.jl"))
include(joinpath("interface","simulate.jl"))
export Network, Simulation, Explicit, Uniform, Mean, Histogram,
  simulate,
  add_species!,   rmv_species!,
  add_reaction!,  rmv_reaction!,
  add_reactant!,  rmv_reactant!,
  add_product!,   rmv_product!,
  add_parameter!, set_parameter!, rmv_parameter!

# Kinetics
include(joinpath("kinetics","mass_action.jl"))

# Algorithms
include(joinpath("algorithm","ssa.jl"))
include(joinpath("algorithm","odm.jl"))
include(joinpath("algorithm","frm.jl"))
include(joinpath("algorithm","nrm.jl"))
include(joinpath("algorithm","sal.jl"))

# Output
include(joinpath("output", "update.jl"))
#include(joinpath("output","pstate.jl"))
#include(joinpath("output","ptrajectory.jl"))
#include(joinpath("output","simresult.jl"))
#include(joinpath("output","simjob.jl"))

end # module
