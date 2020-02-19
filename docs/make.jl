using Documenter, BioSimulator

# setup GR backend for Travis CI
# ENV["GKSwstype"] = "100"
# ENV["PLOTS_TEST"] = "true"

makedocs(
  doctest  = false,
  format   = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
  modules  = [BioSimulator],
  clean    = true,
  sitename = "BioSimulator.jl",
  authors  = "Alfonso Landeros, Timothy Stutz, Kevin L. Keys, Mary E. Sehl, Alexander Alekseyenko, Janet S. Sinsheimer, Kenneth Lange",
  pages = [
    "Home" => "index.md",
    "Model Specification" => "man/model-specification.md",
    "Simulation" => "man/simulation.md",
    "Plotting" => "man/plotting.md",
    "Tutorials" => "man/tutorials.md",
    "Algorithms and References" => "man/algorithms.md",
    "Citing" => "man/citing.md"
  ]
)

deploydocs(
  repo   = "github.com/alanderos91/BioSimulator.jl.git",
)
