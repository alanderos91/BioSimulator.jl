using BioSimulator
using Base.Test

# Load test models
include("test_models.jl")

# List of tests
tests = [#"mass_action",
         #"time_derivatives",
         "network",
         "kendall",
         "linear",
         "independent"]

println("Running tests:")

for t in tests
  println(" * $(t)")
  include("$(t).jl")
end
