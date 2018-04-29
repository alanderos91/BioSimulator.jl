"""
```
AbstractReactionSystem
```

An `AbstractReactionSystem` represents the underlying reactions in a system of interacting particles.

`stoichiometry`: A net stoichiometric matrix that describes how each species (rows) changes when a reaction fires (columns).
`coefficients`: A stoichiometric matrix for the reactants.
`scaled_rates`: A vector of reaction rates, scaled by any additional coefficients (e.g. if the reaction `X + X --> X2` has rate `k`, then its mass-action propensity is `0.5 * k * X * (X - 1)` and the scaled rate is `0.5 * k`.)
`propensities`: A `PropensityVector` that caches each reaction propensity as well as their running sum.
`dependencies`: An adjacency matrix that represents how a given reaction's propensity depends on the other reactions.
"""
abstract type AbstractReactionSystem end

stoichiometry(r::AbstractReactionSystem) = r.stoichiometry
coefficients(r::AbstractReactionSystem)  = r.coefficients
propensities(r::AbstractReactionSystem)  = r.propensities
scaled_rates(r::AbstractReactionSystem)  = r.scaled_rates
dependencies(r::AbstractReactionSystem)  = r.dependencies

@inbounds function update_propensity!{T}(
  a  :: PVec{T},
  r  :: AbstractReactionSystem,
  Xt :: Vector{Int},
  j  :: Integer)

  aⱼ = compute_mass_action(Xt, r, j)
  update_errorbound!(a, aⱼ, j)
  update_intensity!(a,  aⱼ, j)
  a[j] = aⱼ

  return nothing
end

@inbounds function update_all_propensities!{T}(
  a  :: PVec{T},
  r  :: AbstractReactionSystem,
  Xt :: Vector{Int})

  total_sum = zero(T)

  for j in eachindex(a)
    a[j] = compute_mass_action(Xt, r, j)
    total_sum = total_sum + a[j]
  end

  a.intensity   = total_sum
  a.error_bound = length(a) * intensity(a) * eps(T)

  return nothing
end

@inbounds function update_dependent_propensities!{T}(
  a  :: PVec{T},
  r  :: AbstractReactionSystem,
  Xt :: Vector{Int},
  μ  :: Integer)

  dg = dependencies(r)
  dependents = dg[μ]

  for j in dependents
    update_propensity!(a, r, Xt, j)
  end

  return nothing
end

@inline function update_propensities!(
  a  :: PVec,
  r  :: AbstractReactionSystem,
  Xt :: Vector{Int},
  μ  :: Integer)

  if isstable(a)
    update_dependent_propensities!(a, r, Xt, μ)
  else
    update_all_propensities!(a, r, Xt)
  end

  return nothing
end

@inbounds function compute_mass_action(
    Xt :: Vector{Int},
    r  :: AbstractReactionSystem,
    j  :: Integer)

    U = coefficients(r)
    c = scaled_rates(r)

    return c[j] * compute_mass_action(Xt, U, j)
end

@inbounds function compute_mass_action_deriv(
    Xt :: Vector{Int},
    r  :: AbstractReactionSystem,
    j  :: Integer,
    k  :: Integer)

    U = coefficients(r)
    c = scaled_rates(r)

    return c[j] * compute_mass_action_deriv(Xt, U, j, k)
end

fire_reaction!(Xt::Vector{Int}, r::AbstractReactionSystem, μ::Integer) = fire_reaction!(Xt, stoichiometry(r), μ)

function fire_reactions!(Xt::Vector{Int}, r::AbstractReactionSystem, events::Vector{Int})
    for j in eachindex(events)
        fire_reaction!(Xt, stoichiometry(r), j, events[j])
    end
end
