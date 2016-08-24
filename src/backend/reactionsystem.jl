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

  # check for loss of random bits
  if islossy(a)
    fill!(a.cache, zero(eltype(a)))
    a.intensity = zero(eltype(a))
    a.error_bound = zero(eltype(a))

    compute_propensities!(r, Xt)
  else
    dg = dependencies(r)

    dependents = dg[μ]

    for j in dependents
      a[j] = compute_mass_action(Xt, r, j)
    end
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

fire_reaction!(Xt::Vector{Int}, r::AbstractReactionSystem, μ::Integer) = fire_reaction!(Xt, stoichiometry(r), μ)

function fire_reactions!(Xt::Vector{Int}, r::AbstractReactionSystem, events::Vector{Int})
    for j in eachindex(events)
        fire_reaction!(Xt, stoichiometry(r), j, events[j])
    end
end
