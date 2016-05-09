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

# From Cao, Li, & Petzold 2004
function linear(M, x0)
    m = Network("Linear Chain System")

    m <= Species(:S1, x0)

    for i = 2:(M+1)
        m <= Species(symbol("S$(i)"), 0)
    end

    for i = 1:M
        m <= Reaction(symbol("R$(i)"), symbol("k$(i)"),
            r=(symbol("S$(i)") => 1),
            p=(symbol("S$(i+1)") => 1)
        )
        m <= parameter(symbol("k$(i)"), 1.0)
    end

    return m
end

function independent(n, x0)
  m = Network("Independent System")

  for i in 1:n
    m <= Species(symbol("S$(i)"), x0)
    m <= Reaction(symbol("R$(i)"), symbol("k$(i)"), r=(symbol("S$(i)") => 1))
    m <= parameter(symbol("k$(i)"), 1.0)
  end

  return m
end
