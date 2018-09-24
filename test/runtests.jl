using NewBioSimulator
using Test

for model in readdir("models")
  include(joinpath("models", model))
end

include(joinpath("interface", "well-mixed", "network.jl"))

include(joinpath("data-structures", "dep_graph.jl"))
include(joinpath("data-structures", "priority_queue.jl"))
