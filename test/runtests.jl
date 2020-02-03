using BioSimulator
using Test, Statistics

for model in readdir("test-models")
  include(joinpath("test-models", model))
end

include(joinpath("interface", "well-mixed", "network.jl"))

include(joinpath("data-structures", "dep_graph.jl"))
include(joinpath("data-structures", "priority_queue.jl"))

include(joinpath("model", "reaction_system.jl"))

include(joinpath("simulators", "mean_convergence.jl"))

#### spatial simulations
include(joinpath("lattice.jl"))
include(joinpath("ips-interface.jl"))
include(joinpath("execute.jl"))

#### simulation output
include(joinpath("output", "SamplePath.jl"))
