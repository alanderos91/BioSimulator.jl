# Reactions
R1:
  S1 > $pool
  c1 * S1

R2:
  {2} S1 > S2
  c2 * S1 * (S1 - 1) / 2

R3:
  S2 > {2} S1
  c3 * S2

R4:
  S2 > S3
  c4 * S2

# Fixed species

# Variable species
S1 = 1000.0
S2 = 0.0
S3 = 0.0

# Parameters
c1 = 1.0
c2 = 0.002
c3 = 0.5
c4 = 0.04
