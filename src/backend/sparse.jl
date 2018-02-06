"""
```
DenseReactionSystem(reactions, id2ind, c, d)
DenseReactionSystem(stoichiometry, coefficients, scaled_rates, propensities, dependencies)
```

A `SparseReactionSystem` represents the underlying reactions in a system of interacting particles using sparse data structures.

`stoichiometry`: A net stoichiometric matrix that describes how each species (rows) changes when a reaction fires (columns).
`coefficients`: A stoichiometric matrix for the reactants.
`scaled_rates`: A vector of reaction rates, scaled by any additional coefficients (e.g. if the reaction `X + X --> X2` has rate `k`, then its mass-action propensity is `0.5 * k * X * (X - 1)` and the scaled rate is `0.5 * k`.)
`propensities`: A `PropensityVector` that caches each reaction propensity as well as their running sum.
`dependencies`: An adjacency matrix that represents how a given reaction's propensity depends on the other reactions.
"""
struct SparseReactionSystem <: AbstractReactionSystem
  stoichiometry :: SparseMatrixCSC{Int, Int}
  coefficients  :: SparseMatrixCSC{Int, Int}
  scaled_rates  :: Vector{Float64}
  propensities  :: PVec{Float64}
  dependencies  :: Vector{Vector{Int}}
end

function SparseReactionSystem(reactions, id2ind, c, d)
  V, U = sparse_stoichiometry(reactions, id2ind, c, d)

  k = zeros(Float64, d)
  compute_scaled_rates!(k, reactions)

  a = PVec{Float64}(d)

  dg = make_dependency_graph(reactions)

  return SparseReactionSystem(V, U, k, a, dg)
end

function sparse_stoichiometry(reactions, id2ind, c, d)
  V = spzeros(c, d)
  U = spzeros(c, d)

  V, U = compute_net_stoichiometry!(V, U, reactions, id2ind)

  return V, U
end

@inbounds @fastmath function fire_reaction!(
  Xt :: Vector{Int},
  V  :: SparseMatrixCSC{Int, Int},
  μ  :: Integer)

  Vμ  = nonzeros(V)
  idx = rowvals(V)

  for k in nzrange(V, μ)
    i = idx[k]
    Xt[i] = Xt[i] + Vμ[k]
  end

  return nothing
end

@inbounds @fastmath function fire_reaction!(
  Xt :: Vector{Int},
  V  :: SparseMatrixCSC{Int, Int},
  μ  :: Integer,
  n  :: Integer)

  Vμ  = nonzeros(V)
  idx = rowvals(V)

  for k in nzrange(V, μ)
    i = idx[k]
    Xt[i] = Xt[i] + n * Vμ[k]
  end

  return nothing
end

@inbounds @fastmath function compute_mass_action(
  Xt :: Vector{Int},
  U  :: SparseMatrixCSC{Int,Int},
  j  :: Integer)

  rv = rowvals(U)
  nz = nonzeros(U)
  value = 1

  for i in nzrange(U, j)
    value = value * Xt[rv[i]]
    for k in 2:nz[i]
      value = value * (Xt[rv[i]] - (k - 1))
    end
  end

  return value
end

@inbounds @fastmath function compute_mass_action_deriv(
  Xt :: Vector{Int},
  U  :: SparseMatrixCSC{Int,Int},
  j  :: Integer,
  k  :: Integer)

  rv = rowvals(U)
  nz = nonzeros(U)
  value = 0

  if U[k,j] > 0
    value = 1
    for i in nzrange(U, j)
      if rv[i] != k
        value = value * Xt[rv[i]]
        for n in 2:nz[i]
          value = value * (Xt[rv[i]] - (n - 1))
        end
      else
        tmp1 = 0
        for m in 1:nz[i]
          tmp2 = 1
          for n in 1:nz[i]
            if m != n; tmp2 = tmp2 * (Xt[rv[i]] - (n - 1)); end
          end
          tmp1 = tmp1 + tmp2
        end
        value = value * tmp1
      end
    end
  end

  return value
end
