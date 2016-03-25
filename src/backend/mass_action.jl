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
