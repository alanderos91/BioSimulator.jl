using NewBioSimulator
using Test

for model in readdir("models")
  include(joinpath("models", model))
end

include("network.jl")
include("dep_graph.jl")
