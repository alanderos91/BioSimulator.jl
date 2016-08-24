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

@inbounds function update_propensity!(r::AbstractReactionSystem, Xt, j::Integer)
  a = propensities(r)

  a_j = compute_mass_action(Xt, r, j)
  update_errorbound!(a, a_j, j)
  update_intensity!(a,  a_j, j)
  a[j] = a_j

  return r
end

@inbounds function compute_propensities!(r::AbstractReactionSystem, Xt, iterable::AbstractVector)

  for j in iterable
    update_propensity!(r, Xt, j)
  end

  return r
end

function update_all_propensities!(r::AbstractReactionSystem, Xt)
  a = propensities(r)

  compute_propensities!(r, Xt, eachindex(a))

  return r
end

function update_dependent_propensities!(r::AbstractReactionSystem, Xt, μ::Integer)
  a  = propensities(r)
  dg = dependencies(r)

  compute_propensities!(r, Xt, dg[μ])

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
