# Reactions
receptor_upregulation:
  $pool > R
  k1

receptor_downregulation:
  R > $pool
  k2 * R

ligand_binding:
  L + R > RL + L
  k3 * L * R

ligand_degradation:
  RL > R
  k4 * RL

protein_activation:
  RL + G > Ga + Gbg + RL
  k5 * RL * G

dephosphorylation:
  Ga > Gd
  k6 * Ga

rebinding:
  Gd + Gbg > G
  k7 * Gd * Gbg

bound_receptor_upregulation:
  $pool > RL
  k8

# Fixed species

# Variable species
R = 50
L = 2
RL = 0
G = 50
Ga = 0
Gbg = 0
Gd = 0

# Parameters
k1 = 0.038
k2 = 0.0004
k3 = 0.042
k4 = 0.01
k5 = 0.011
k6 = 0.1
k7 = 0.00105
k8 = 3.21
