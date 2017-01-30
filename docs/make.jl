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
    "Interface"  => "interface.md",
    "Algorithms" => "algorithms.md"
  ]
)

deploydocs(
  repo   = "bitbucket.com/alanderos/biosimulator.jl.git",
  target = "build",
  deps   = nothing,
  make   = nothing
)
