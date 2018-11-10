function yeast(;
  R=50, L=2, RL=0, G=50, Ga=0, Gbg=0, Gd=0,
  k1=0.0038, k2=4e-4, k3=0.042, k4=0.01, k5=0.011, k6=0.1, k7=1.05e-3, k8=3.21)
  m = Network("Yeast polarization")

  m <= Species("R", R) # receptor
  m <= Species("L", L) # ligand
  m <= Species("RL", RL) # receptor-ligand
  m <= Species("G", G) # G protein
  m <= Species("Ga", Ga) # G alpha subunit
  m <= Species("Gbg", Gbg) # G beta-gamma subunit
  m <= Species("Gd", Gd) # G delta subunit

  m <= Reaction("receptor upregulation", k1, "0 --> R")
  m <= Reaction("receptor downregulation", k2, "R --> 0")
  m <= Reaction("ligand binding", k3, "L + R --> RL + L")
  m <= Reaction("ligand degradation", k4, "RL --> R")
  m <= Reaction("protein activation", k5, "RL + G --> Ga + Gbg + RL")
  m <= Reaction("dephosphorylation", k6, "Ga --> Gd")
  m <= Reaction("rebinding", k7, "Gd + Gbg --> G")
  m <= Reaction("bound receptor upregulation", k8, "0 --> RL")
  
  return m
end

model = yeast()
