export Reaction

type Reaction
  id::ASCIIString
  rate::ASCIIString
  propensity::Float64
  pre::Vector{Int}
  post::Vector{Int}

  function Reaction(id::ASCIIString, rate::ASCIIString, pre::Vector{Int}, post::Vector{Int})
    if any(pre .< 0) || any(post .< 0)
      error("Stoichiometric coefficients must be positive.")
    end
    return new(id, rate, 0.0, pre, post)
  end
end
