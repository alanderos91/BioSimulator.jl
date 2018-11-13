function F_mmek(x, parms)
  S, E, SE, P = x
  k1, k2, k3 = parms

  binding = k1 * S * E
  dissociation = k2 * SE
  conversion = k3 * SE

  return [binding, dissociation, conversion]
end

function gillespie_mmek(;
  S = 301, E = 120, SE = 0, P = 0,
  k1 = 0.00166, k2 = 0.0001, k3 = 0.1)
  # initial conditions
  x0 = [S, E, SE, P]
  # stoichiometries
  nu = [
    [-1 -1 1 0]; # binding, S + E --> SE
    [1 1 -1 0];  # dissociation, SE --> S + E
    [0 1 -1 1]   # conversion, SE --> P + E
  ]
  # parameter vector
  parms = [k1, k2, k3]

  return x0, nu, parms
end

x0, nu, parms = gillespie_mmek()
F = F_mmek
