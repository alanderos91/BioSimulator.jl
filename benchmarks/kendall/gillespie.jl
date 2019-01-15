function F_kendall(x, parms)
  X = x[1]
  α, μ, ν = parms

  birth = α * X
  death = μ * X
  immigration = ν

  return [birth, death, immigration]
end

function gillespie_kendall(;
  X = 5,
  α = 2.0, μ = 1.0, ν = 0.5)
  # initial conditions
  x0 = [X]
  # stoichiometries:
  #   birth, X --> X + X 
  #   death, X --> 0
  #   immigration, 0 --> X
  nu = [1 -1 1]
  # parameter vector
  parms = [α, μ, ν]

  return x0, nu, parms
end

x0, nu, parms = gillespie_kendall()
F = F_kendall
