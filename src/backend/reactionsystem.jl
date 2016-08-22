##### AbstractReactionSystem #####
#=
  stoichiometry :: T
  coefficients  :: T
  scaled_rates  :: Vector{Float64}
  propensities  :: PVec{Float64}
  dependencies  :: Vector{Vector{Int}}
=#

abstract AbstractReactionSystem

stoichiometry(r::AbstractReactionSystem) = r.stoichiometry
coefficients(r::AbstractReactionSystem)  = r.coefficients
propensities(r::AbstractReactionSystem)  = r.propensities
scaled_rates(r::AbstractReactionSystem)  = r.scaled_rates
dependencies(r::AbstractReactionSystem)  = r.dependencies

function compute_propensities!(r::AbstractReactionSystem, Xt)
  a = propensities(r)

  for j in eachindex(a)
    a[j] = compute_mass_action(Xt, r, j)
  end

  return r
end

function update_propensities!(r::AbstractReactionSystem, Xt, μ)
  a  = propensities(r)
  dg = dependencies(r)

  dependents = dg[μ]

  for j in dependents
    a[j] = compute_mass_action(Xt, r, j)
  end

  return r
end

function compute_mass_action(
    Xt :: Vector{Int},
    r  :: AbstractReactionSystem,
    j  :: Integer)

    U = coefficients(r)
    c = scaled_rates(r)

    return c[j] * compute_mass_action(Xt, U, j)
end

function compute_mass_action_deriv(
    Xt :: Vector{Int},
    r  :: AbstractReactionSystem,
    j  :: Integer,
    k  :: Integer)

    U = coefficients(r)
    c = scaled_rates(r)

    return c[j] * compute_mass_action_deriv(Xt, U, j, k)
end

fire_reaction!(Xt, r, μ) = fire_reaction!(Xt, stoichiometry(r), μ)

function fire_reactions!(Xt, r, events)
    for j in eachindex(events)
        fire_reaction!(Xt, stoichiometry(r), j, events[j])
    end
end
