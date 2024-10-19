mutable struct PoissonLeapMethod{F1,F2} <: UnsafeLeapAlgorithm
  rates::Vector{Float64}
  total_rate::Float64
  stoichiometry::SparseMatrixCSC{Int,Int}
  leap_formula::F1
  validate_leap!::F2
end

function generate_events!(::PoissonLeapMethod, jumps, rates, s)
  if s < 0
    throw(ArgumentError("The leap length $(s) should be positive; got s = $(s)."))
  end
  map!(rate -> pois_rand(rate * s), jumps, rates)
  return nothing
end

function update!(algorithm::PoissonLeapMethod, state, model)
  # unpack
  rates = algorithm.rates
  total_rate = algorithm.total_rate
  stoichiometry = algorithm.stoichiometry

  # update jump rates
  update_jump_rates!(algorithm, state, model)

  # update auxiliary variables
  update!(algorithm.leap_formula, state, model, stoichiometry, rates, total_rate)

  return nothing
end
