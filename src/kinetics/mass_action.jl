function mass_action_helper(stoich::Vector{Int}, x::Vector{Species})
  acc = 1.0

  for i in eachindex(stoich)
    for j = 1:stoich[i]
      acc = acc * (x[i].pop - (j - 1))
    end
  end

  return acc
end

function mass_action_helper(stoich::Vector{Int}, x::Vector{Species}, k::Int)
  acc = 1.0
  for i in eachindex(stoich)
    if i != k
      for j = 1:stoich[i]
        acc = acc * (x[i].pop - (j - 1))
      end
    end
  end
  return acc
end

function mass_action(c::Float64, stoich::Vector{Int}, x::Vector{Species})
  return c * mass_action_helper(stoich, x)
end

function mass_action(r::Reaction, x::Vector{Species}, params::Dict{ASCIIString, Float64})
  c = params[r.rate]
  return mass_action(c, r.pre, x)
end

function mass_action_deriv(c::Float64, stoich::Vector{Int}, x::Vector{Species}, k::Int)
  acc = 1.0
  if stoich[k] == 0
    acc = 0.0
  else
    acc = mass_action_helper(stoich, x, k)
  end

  if stoich[k] == 2
    acc = acc * (2 * x[k].pop - 1)
  elseif stoich[k] > 2
    error("Higher order reactions are not supported.")
  end

  return c * acc
end

function mass_action_deriv(r::Reaction, x::Vector{Species}, params::Dict{ASCIIString, Float64}, k::Int)
  c = params[r.rate]
  return mass_action_deriv(c, r.pre, x, k)
end

function propensity!(r::Reaction, x::Vector{Species}, param::Dict{ASCIIString,Float64})
  r.propensity = mass_action(r, x, param)
  return;
end
