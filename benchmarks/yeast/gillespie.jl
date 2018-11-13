function F_yeast(x, parms)
  R, L, RL, G, Ga, Gbg, Gd = x
  k1, k2, k3, k4, k5, k6, k7, k8 = parms

  receptor_uregulate = k1
  receptor_dregulate = k2 * R
  ligand_binding = k3 * L * R
  ligand_degrade = k4 * RL
  G_activation = k5 * RL * G
  G_dephosphorylate = k6 * Ga
  G_rebinding = k7 * Gd * Gbg
  bound_receptor_uregulate = k8

  rates = [
    receptor_uregulate,
    receptor_dregulate,
    ligand_binding,
    ligand_degrade,
    G_activation,
    G_dephosphorylate,
    G_rebinding,
    bound_receptor_uregulate
  ]
  return rates
end

function gillespie_yeast(;
  R=50, L=2, RL=0, G=50, Ga=0, Gbg=0, Gd=0,
  k1=0.0038, k2=4e-4, k3=0.042, k4=0.01, k5=0.011, k6=0.1, k7=1.05e-3, k8=3.21)
  # initial conditions
  x0 = [R, L, RL, G, Ga, Gbg, Gd]
  # stoichiometries
  nu = [
    [1 0 0 0 0 0 0];
    [-1 0 0 0 0 0 0];
    [-1 0 1 0 0 0 0];
    [1 0 -1 0 0 0 0];
    [0 0 0 -1 1 1 0];
    [0 0 0 0 -1 1 0];
    [0 0 0 1 0 -1 -1];
    [0 0 1 0 0 0 0]
  ]
  # parameter vector
  parms = [k1, k2, k3, k4, k5, k6, k7, k8]

  return x0, nu, parms
end

x0, nu, parms = gillespie_yeast()
F = F_yeast