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

function autoreg()
  m = Network("auto-regulation")

  m <= Species("gene",   10)
  m <= Species("P2_gene", 0)
  m <= Species("RNA",     0)
  m <= Species("P",       0)
  m <= Species("P2",      0)

  m <= Reaction("repression binding",         :k1,  :(gene + P2 --> P2_gene))
  m <= Reaction("reverse repression binding", :k1r, :(P2_gene --> gene + P2))
  m <= Reaction("transcription",              :k2,  :(gene --> gene + RNA))
  m <= Reaction("translation",                :k3,  :(RNA --> RNA + P))
  m <= Reaction("dimerization",               :k4,  :(P + P --> P2))
  m <= Reaction("dissociation",               :k4r, :(P2 --> P + P))
  m <= Reaction("RNA degradation",            :k5,  :(RNA --> 0))
  m <= Reaction("protein degradation",        :k6,  :(P --> 0))

  m <= Parameter(:k1,   1.0)
  m <= Parameter(:k1r, 10.0)
  m <= Parameter(:k2,  0.01)
  m <= Parameter(:k3,  10.0)
  m <= Parameter(:k4,   1.0)
  m <= Parameter(:k4r,  1.0)
  m <= Parameter(:k5,   0.1)
  m <= Parameter(:k6,  0.01)

  return m
end
