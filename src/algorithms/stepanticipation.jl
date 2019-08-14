mutable struct StepAnticipationMethod{F1,F2} <: UnsafeLeapAlgorithm
  rates::Vector{Float64}
  total_rate::Float64
  leap_formula::F1
  validate_leap!::F2
  U::SparseMatrixCSC{Int,Int} # reactant stoichiometry
  V::SparseMatrixCSC{Int,Int} # net stoichiometry
  dxdt::Vector{Float64}
  drdt::Vector{Float64}
end

function generate_events!(f::StepAnticipationMethod, v, rates, s)
  map!((rate, drdt) -> pois_rand(rate * s + 1//2 * drdt * s * s), v, rates, f.drdt)
end

function update!(algorithm::StepAnticipationMethod, state, model)
  # unpack
  rates = algorithm.rates
  total_rate = algorithm.total_rate
  V = algorithm.V
  U = algorithm.U

  dxdt = algorithm.dxdt
  drdt = algorithm.drdt

  # update jump rates
  update_jump_rates!(algorithm, state, model)

  # update auxiliary variables
  update!(algorithm.leap_formula, state, model, V, rates, total_rate)

  return nothing
end

