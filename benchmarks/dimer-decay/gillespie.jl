function F_dimerdecay(x, parms)
  S1, S2, S3 = x
  c1, c2, c3, c4 = parms

  R1 = c1 * S1
  R2 = c2 * S1 * (S1 - 1) / 2
  R3 = c3 * S2
  R4 = c4 * S2

  return [R1, R2, R3, R4]
end

function gillespie_dimerdecay(;
  x = [1_000, 0, 0], c = [1.0, 0.002, 0.5, 0.04]
  )
  # initial conditions
  x0 = copy(x)
  # stoichiometries
  nu = [
    [-1  0  0]; # S1 --> 0
    [-2  1  0]; # 2 * S1 --> S2
    [ 2 -1  0]; # S2 --> 2 * S1
    [ 0 -1  1]  # S2 --> S3
  ]
  # parameter vector
  parms = copy(c)
  return x0, nu, parms
end

x0, nu, parms = gillespie_dimerdecay()
F = F_dimerdecay
