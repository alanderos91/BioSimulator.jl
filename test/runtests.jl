using BioSimulator
using Base.Test

import BioSimulator: Algorithm

function run_test{T}(
  model :: Network,
  alg   :: T,
  t     :: Real,
  n     :: Integer,
  m     :: Integer;
  kwargs...)
  result = simulate(model, alg, time=t, epochs=n, trials=m, kwargs...)
end

# Load test models
include("test_models.jl")

# List of tests
tests = ["mass_action",
         "time_derivatives",
         "network",
         "sort",
         #"kendall",
         "linear",
         #"independent",
         "sir"]

println("Running tests:")

for t in tests
  println(" * $(t)")
  include("$(t).jl")
end
