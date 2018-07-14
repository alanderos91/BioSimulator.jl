using Documenter, BioSimulator

makedocs(
  doctest  = false,
  format   = :html,
  modules  = [BioSimulator],
  clean    = true,
  sitename = "BioSimulator.jl",
  authors  = "Alfonso Landeros, Mary E. Sehl, Kenneth Lange",
  pages = [
    "Home"       => "index.md",
    "Overview"   => "man/overview.md",
    "Modeling"   => "man/modeling.md",
    "Algorithms" => "man/algorithms.md",
    "Examples"   => "man/examples.md",
  ]
)

deploydocs(
  repo   = "github.com/alanderos91/BioSimulator.jl.git",
  target = "build",
  julia  = "0.6",
  deps   = nothing,
  make   = nothing
)
