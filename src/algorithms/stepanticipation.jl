struct StepAnticipationMethod{F1,F2} <: UnsafeLeapAlgorithm
  rates::Vector{Float64}
  total_rate::Float64
  stoichiometry::SparseMatrixCSC{Int,Int}
  dxdt::Vector{Float64}
  drdt::Vector{Float64}
  leap_formula::F1
  validate_leap!::F2
end

# function initialize!(algorithm::StepAnticipationMethod, state, model)
#   update!(algorithm, state, model)
# end

function generate_events!(f::StepAnticipationMethod, v, rates, s)
  map!((rate, drdt) -> pois_rand(rate * s + 1//2 * drdt * s * s), v, rates, f.drdt)
end

function update!(algorithm::StepAnticipationMethod, state, model)
  # unpack
  rates = algorithm.rates
  total_rate = algorithm.total_rate
  stoichiometry = algorithm.stoichiometry

  # update jump rates
  update_jump_rates!(algorithm, state, model)

  # update auxiliary variables
  update!(algorithm.leap_formula, state, model, stoichiometry, rates, total_rate)
  mean_derivatives!()
  time_derivatives!()
end

function mean_derivatives!(dxdt, rates, V)
  nz = nonzeros(V)
  ix = rowvals(V)

  fill!(dxdt, 0.0)

  for j in eachindex(rates)
    for k in nzrange(V, j)
      i = ix[k]
      dxdt[i] = dxdt[i] + rates[j] * nz[k]
    end
  end

  return dxdt
end

function time_derivatives!(drdt, x, U, dxdt)
  # U - reactant coefficients!
  ix = rowvals(U)

  for j in eachindex(drdt)
    drdt[j]  = 0.0
    for k in nzrange(U, j)
      drdx = rate_derivative(model, state, ix[k], j)
      drdt[j] = drdt[j] + drdx * dxdt[ix[k]]
    end
  end

  return drdt
end
