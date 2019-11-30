function kendall(; x0 = 5, birth_rate = 2.0, death_rate = 1.0, immgr_rate = 0.5)
  m = Network("Kendall's Process")

  m <= Species("X", x0)

  m <= Reaction("birth",       birth_rate, "X --> X + X")
  m <= Reaction("death",       death_rate, "X --> 0")
  m <= Reaction("immigration", immgr_rate, "0 --> X")
  
  return m
end
