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
    "Benchmarks" => [],
    "Developers" => []
  ]
)

deploydocs(
  repo   = "bitbucket.com/alanderos/biosimulator.jl.git",
  target = "build",
  deps   = Deps.pip("mkdocs", "python-markdown-math"),
  make   = nothing
)
