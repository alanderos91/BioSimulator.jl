using Documenter, BioSimulator

# setup GR backend for Travis CI
ENV["GKSwstype"] = "100"
ENV["PLOTS_TEST"] = "true"

makedocs(
  doctest  = false,
  format   = Documenter.HTML(),
  modules  = [BioSimulator],
  clean    = true,
  sitename = "BioSimulator.jl",
  authors  = "Alfonso Landeros, Timothy Stutz, Kevin L. Keys, Mary E. Sehl, Alexander Alekseyenko, Janet S. Sinsheimer, Kenneth Lange",
  pages = [
    "Home"       => "index.md",
    "Algorithms" => "man/algorithms.md",
  ]
)

deploydocs(
  repo   = "github.com/alanderos91/BioSimulator.jl.git",
  target = "build",
  deps   = nothing,
  make   = nothing
)
