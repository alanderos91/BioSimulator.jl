function kendall()
  m = Network("Kendall's Process")

  m <= Species(:X, 5)

  m <= Reaction(:Birth,       :α, :(X --> X + X))
  m <= Reaction(:Death,       :μ, :(X --> 0))
  m <= Reaction(:Immigration, :ν, :(0 --> X))

  m <= Parameter(:α, 2.0)
  m <= Parameter(:μ, 1.0)
  m <= Parameter(:ν, 0.5)

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
        m <= Parameter(symbol(:k, i), 1.0)
    end

    return m
end

function independent(n, x0)
  m = Network("Independent System")

  for i in 1:n
    m <= Species(symbol(:S, i), x0)
    m <= Reaction(symbol(:R, i), symbol(:k, i), Expr(:-->, symbol(:S,i), 0))
    m <= Parameter(symbol(:k, i), 1.0)
  end

  return m
end
