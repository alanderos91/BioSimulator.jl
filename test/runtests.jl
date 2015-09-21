using BioSimulator
using Base.Test

tests = ["mass_action",
         "species",
         "time_derivatives"]

println("Running tests:")

for t in tests
  println(" * $(t)")
  include("$(t).jl")
end
