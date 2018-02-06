"""
```
DenseReactionSystem(reactions, id2ind, c, d)
DenseReactionSystem(stoichiometry, coefficients, scaled_rates, propensities, dependencies)
```

A `DenseReactionSystem` represents the underlying reactions in a system of interacting particles using dense data structures.

`stoichiometry`: A net stoichiometric matrix that describes how each species (rows) changes when a reaction fires (columns).
`coefficients`: A stoichiometric matrix for the reactants.
`scaled_rates`: A vector of reaction rates, scaled by any additional coefficients (e.g. if the reaction `X + X --> X2` has rate `k`, then its mass-action propensity is `0.5 * k * X * (X - 1)` and the scaled rate is `0.5 * k`.)
`propensities`: A `PropensityVector` that caches each reaction propensity as well as their running sum.
`dependencies`: An adjacency matrix that represents how a given reaction's propensity depends on the other reactions.
"""
struct DenseReactionSystem <: AbstractReactionSystem
  stoichiometry :: Matrix{Int}
  coefficients  :: Matrix{Int}
  scaled_rates  :: Vector{Float64}
  propensities  :: PVec{Float64}
  dependencies  :: Vector{Vector{Int}}
end

function DenseReactionSystem(reactions, id2ind, c, d)
  V, U = dense_stoichiometry(reactions, id2ind, c, d)

  k = zeros(Float64, d)
  compute_scaled_rates!(k, reactions)

  a = PVec{Float64}(d)

  dg = make_dependency_graph(reactions)

  return DenseReactionSystem(V, U, k, a, dg)
end

function dense_stoichiometry(reactions, id2ind, c, d)
  V = zeros(c, d)
  U = zeros(c, d)

  V, U = compute_net_stoichiometry!(V, U, reactions, id2ind)

  return V, U
end

@inbounds @fastmath function fire_reaction!(
  Xt :: Vector{Int},
  V  :: Matrix{Int},
  μ  :: Integer)

  for k in eachindex(Xt)
    Xt[k] = Xt[k] + V[k, μ]
  end

  return nothing
end

@inbounds @fastmath function fire_reaction!(
  Xt :: Vector{Int},
  V  :: Matrix{Int},
  μ  :: Integer,
  n  :: Integer)

  for k in eachindex(Xt)
    Xt[k] = Xt[k] + n * V[k, μ]
  end

  return nothing
end

@inbounds @fastmath function compute_mass_action(
  Xt :: Vector{Int},
  U  :: Matrix{Int},
  j  :: Integer)

  value = 1

  for i in eachindex(Xt)
    if U[i, j] > 0
      value = value * Xt[i]
      for k in 2:U[i, j]
        value = value * (Xt[i] - (k - 1) )
      end
    end
  end

  return value
end

@inbounds @fastmath function compute_mass_action_deriv(
  Xt :: Vector{Int},
  U  :: Matrix{Int},
  j  :: Integer,
  k  :: Integer)

  value = 0

  if U[k, j] > 0
    value = 1
    for i in eachindex(Xt)
      if i != k
        if U[i, j] > 0
          value = value * Xt[i]
          for n in 2:U[i, j]
            value = value * (Xt[i] - (n - 1))
          end
        end
      elseif U[i, j] > 1
        tmp1 = 0
        for m in 1:U[i, j]
          tmp2 = 1
          for n in 1:U[i, j]
            if m != n; tmp2 = tmp2 * (Xt[i] - (n - 1)); end
          end
          tmp1 = tmp1 + tmp2
        end
        value = value * tmp1
      end
    end
  end

  return value
end
