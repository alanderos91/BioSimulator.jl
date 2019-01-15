function mmek(;
  S=301, E=120, SE=0, P=0,
  k1=0.00166, k2=0.0001, k3=0.1)
  # Initialize network
  m = Network("Michaelis-Menten")
  
  # Species Definitions
  m <= Species("S",  S)
  m <= Species("E",  E)
  m <= Species("SE", SE)
  m <= Species("P",  P)
  
  # Reaction Definitions
  m <= Reaction("binding",      k1, "S + E --> SE")
  m <= Reaction("dissociation", k2, "SE --> S + E")
  m <= Reaction("conversion",   k3, "SE --> E + P")
  
  return m
end

model = mmek()
