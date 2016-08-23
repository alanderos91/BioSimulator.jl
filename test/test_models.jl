function kendall(X0, α, μ, ν)
  m = Network("Kendall's Process")

  m <= Species(:X, X0)

  m <= Reaction(:Birth,       α, :(X --> X + X))
  m <= Reaction(:Death,       μ, :(X --> 0))
  m <= Reaction(:Immigration, ν, :(0 --> X))

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
        m <= Reaction(symbol(:R, i), 1.0, Expr(:-->, symbol(:S,i), symbol(:S,i+1)))
    end

    return m
end

function independent(n, x0)
  m = Network("Independent System")

  for i in 1:n
    m <= Species(symbol(:S, i), x0)
    m <= Reaction(symbol(:R, i), 1.0, Expr(:-->, symbol(:S,i), 0))
  end

  return m
end

function autoreg(;k1=1.0, k1r=10.0, k2=0.01, k3=10.0, k4=1.0, k4r=1.0, k5=0.1, k6=0.01)
  m = Network("auto-regulation")

  m <= Species("gene",   10)
  m <= Species("P2_gene", 0)
  m <= Species("RNA",     0)
  m <= Species("P",       0)
  m <= Species("P2",      0)

  m <= Reaction("repression binding",         k1,  :(gene + P2 --> P2_gene))
  m <= Reaction("reverse repression binding", k1r, :(P2_gene --> gene + P2))
  m <= Reaction("transcription",              k2,  :(gene --> gene + RNA))
  m <= Reaction("translation",                k3,  :(RNA --> RNA + P))
  m <= Reaction("dimerization",               k4,  :(P + P --> P2))
  m <= Reaction("dissociation",               k4r, :(P2 --> P + P))
  m <= Reaction("RNA degradation",            k5,  :(RNA --> 0))
  m <= Reaction("protein degradation",        k6,  :(P --> 0))

  return m
end
