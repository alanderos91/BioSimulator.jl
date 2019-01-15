function dimerdecay(;x = [1_000, 0, 0], c = [1.0, 0.002, 0.5, 0.04])
  m = Network("dimer-decay / Rathiam, Petzold, Cao, Gillespie 2003")

  m <= Species("S1", x[1])
  m <= Species("S2", x[2])
  m <= Species("S3", x[3])

  m <= Reaction("R1", c[1], "S1 --> 0")
  m <= Reaction("R2", c[2], "2 * S1 --> S2")
  m <= Reaction("R3", c[3], "S2 --> 2 * S1")
  m <= Reaction("R4", c[4], "S2 --> S3")

  return m
end

model = dimerdecay()
