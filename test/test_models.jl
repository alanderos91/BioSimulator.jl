function kendall()
  m = Network("Kendall's Process")

  m <= Species(:X, 5)

  m <= Reaction(:Birth,       :α, :(X --> X + X))
  m <= Reaction(:Death,       :μ, :(X --> 0))
  m <= Reaction(:Immigration, :ν, :(0 --> X))

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

  m <= Reaction(:zero,     :k0,  :(0 --> X))
  m <= Reaction(:first,    :k1,  :(X --> 0))
  m <= Reaction(:second_a, :k2a, :(X + Y --> XY))
  m <= Reaction(:second_b, :k2b, :(X + X --> Y))

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
        m <= Species(symbol(:S, i), 0)
    end

    for i = 1:M
        m <= Reaction(symbol(:R, i), symbol(:k, i), Expr(:-->, symbol(:S,i), symbol(:S,i+1)))
        m <= parameter(symbol(:k, i), 1.0)
    end

    return m
end

function independent(n, x0)
  m = Network("Independent System")

  for i in 1:n
    m <= Species(symbol(:S, i), x0)
    m <= Reaction(symbol(:R, i), symbol(:k, i), Expr(:-->, symbol(:S,i), 0))
    m <= parameter(symbol(:k, i), 1.0)
  end

  return m
end
