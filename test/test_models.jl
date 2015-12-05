function kendall()
  m = Network("Kendall's Process")

  m <= Species(:X, 5)

  m <= Reaction(:Birth,       :α, r=(:X => 1), p=(:X => 2))
  m <= Reaction(:Death,       :μ, r=(:X => 1)             )
  m <= Reaction(:Immigration, :ν,              p=(:X => 1))

  m <= parameter(:α, 2.0)
  m <= parameter(:μ, 1.0)
  m <= parameter(:ν, 0.5)

  return m
end

function test1()
  m = Network("test1")

  m <= Species(:X, 100)
  m <= Species(:Y, 500)
  m <= Species(:XY,  1)

  m <= Reaction(:zero,     :k0,  p=(:X => 1))
  m <= Reaction(:first,    :k1,  r=(:X => 1))
  m <= Reaction(:second_a, :k2a, r=(:X => 1, :Y => 1), p=(:XY => 1))
  m <= Reaction(:second_b, :k2b, r=(:X => 2),          p=(:Y => 1))

  m <= parameter(:k0,  0.1)
  m <= parameter(:k1,  0.01)
  m <= parameter(:k2a, 0.25)
  m <= parameter(:k2b, 0.5)

  return m
end
