function kendall(;
  α :: Float64 = 2.0,
  μ :: Float64 = 1.0,
  ν :: Float64 = 0.5)
  m = Network("Kendall's Process")

  m <= Species("X", 5)

  m <= Reaction("birth", α, "X --> X + X")
  m <= Reaction("death", μ, "X --> 0")
  m <= Reaction("immigration", ν, "0 --> X")

  return m
end

model = kendall()
