using BioSimulator
using Test
using Random
using Statistics

import BioSimulator: Algorithm
import Random: seed!
import Printf: @printf

function run_test(
  model :: Network,
  alg   :: T,
  t     :: Real,
  n     :: Integer,
  m     :: Integer;
  kwargs...) where T
  result = simulate(model, alg, time=t, epochs=n, trials=m, kwargs...)
end

# Load test models
include("test_models.jl")

# List of tests
tests = ["mass_action",
         "time_derivatives",
         "network",
         "sort",
         "summary",
         "kendall",
         "linear",
         "independent",
         "sir"]

println("Running tests:")

for t in tests
  @testset "$(t)" begin
    include("$(t).jl")
  end
end
