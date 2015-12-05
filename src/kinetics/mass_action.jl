function mass_action_helper(stoich::Vector{Int}, x::Vector{Int})
  acc = 1.0

  for i in eachindex(stoich)
    for j = 1:stoich[i]
      acc = acc * (x[i] - (j - 1))
    end
  end

  return acc
end

function mass_action_helper(stoich::Vector{Int}, x::Vector{Int}, k::Int)
  acc = 1.0
  for i in eachindex(stoich)
    if i != k
      for j = 1:stoich[i]
        acc = acc * (x[i] - (j - 1))
      end
    end
  end
  return acc
end

function mass_action(c::Parameter, stoich::Vector{Int}, x::Vector{Int})
  return c * mass_action_helper(stoich, x)
end

function mass_action(r::ReactionChannel, x::Vector{Int}, params::Parameters)
  c = params[r.rate]
  return mass_action(c, r.pre, x)
end

function mass_action_deriv(c::Parameter, stoich::Vector{Int}, x::Vector{Int}, k::Int)
  acc = 1.0
  if stoich[k] == 0
    acc = 0.0
  else
    acc = mass_action_helper(stoich, x, k)
  end

  if stoich[k] == 2
    acc = acc * (2 * x[k] - 1)
  elseif stoich[k] > 2
    error("Higher order reactions are not supported.")
  end

  return c * acc
end

function mass_action_deriv(r::ReactionChannel, x::Vector{Int}, params::Parameters, k::Int)
  c = params[r.rate]
  return mass_action_deriv(c, r.pre, x, k)
end

function propensity!(r::ReactionChannel, x::Vector{Int}, param::Parameters)
  r.propensity = mass_action(r, x, param)
  return;
end

function compute_propensities!(rxns::ReactionVector, x::Vector{Int}, param::Parameters)
  intensity = 0.0
  for r in rxns
    propensity!(r, x, param)
    intensity = intensity + r.propensity
  end
  return intensity
end
