"""
```
OTL
```

Ordinary τ-leaping

"""
mutable struct OTL <: TauLeapMethod
  # parameters
  end_time :: Float64
  ϵ :: Float64
  δ :: Float64
  β :: Float64

  # state variables
  t :: Float64

  # cache
  ξ :: Vector{Float64}
  b :: Matrix{Float64}
  τ :: Vector{Float64}
  events :: Vector{Int}

  # statistics
  stats :: Dict{Symbol,Int}

  function OTL(end_time, ϵ, δ, β)
    new(end_time, ϵ, δ, β,
    0.0, Float64[], Matrix{Float64}(0,0), Float64[], Int[],
    Dict{Symbol,Int}(
      :negative_excursions => 0,
      :contractions => 0,
      :leaping_steps => 0,
      :gillespie_steps => 0
    ))
  end
end

OTL(end_time; ϵ::Float64 = 0.125, δ::Float64 = 100.0, β::Float64 = 0.75, na...) = OTL(end_time, ϵ, δ, β)

set_time!(algorithm::OTL, τ) = (algorithm.t = algorithm.t + τ)

function init!(x :: OTL, Xt, r)
  c = length(Xt)
  d = size(stoichiometry(r), 2)

  x.ξ = zeros(c)    # store expected change
  x.b = zeros(c, d) # store partial derivatives
  x.τ = zeros(d)    # compute leap size
  x.events = zeros(Int, d)

  return nothing
end

function reset!(x :: OTL, a :: PVec)
  x.t = zero(x.t)

  return nothing
end

function step!(algorithm :: OTL, Xt, r)
  a = propensities(r)
  iscritical = (intensity(a) < algorithm.δ)

  if intensity(a) > 0
    if iscritical
      # Gillespie update
      algorithm.stats[:gillespie_steps] += 1
      τ = randexp() / intensity(a)
      set_time!(algorithm, τ)

      if !done(algorithm)
        μ = select_reaction(a)
        fire_reaction!(Xt, r, μ)
        update_propensities!(a, r, Xt, μ)
      end
    else
      # τ-leap update
      algorithm.stats[:leaping_steps] += 1
      expected_change!(algorithm.ξ, a.cache, r)
      propensity_derivatives!(algorithm.b, Xt, r)

      # compute leap size and check for end of interval
      τ = tauleap_DG2001(algorithm.ξ, algorithm.b, algorithm.τ, a.intensity, algorithm.ϵ)
      τ = min(τ, end_time(algorithm) - get_time(algorithm))

      # compute jumps
      for j in eachindex(algorithm.events)
        algorithm.events[j] = poisrand(τ * a[j])
      end

      isbadleap = is_badleap(Xt, r, algorithm.events)
      if isbadleap
        algorithm.stats[:negative_excursions] += 1
      end
    
      while isbadleap
        algorithm.stats[:contractions] += 1
        contract!(algorithm.events, algorithm.β)
        τ = τ * algorithm.β
        isbadleap = is_badleap(Xt, r, algorithm.events)
      end
    
      fire_reactions!(Xt, r, algorithm.events)

      # state dependent updates; these could be fused together
      update_all_propensities!(a, r, Xt)
      # expected_change!(algorithm.ξ, a.cache, r)
      # propensity_derivatives!(algorithm.b, Xt, r)

      set_time!(algorithm, τ)
    end
  elseif intensity(a) == 0
    algorithm.t = algorithm.end_time
  else
    throw(error("intensity = $(intensity(a)) < 0 at time $algorithm.t"))
  end

  return nothing
end

"""
  Expected change per unit time.
"""
function expected_change!(ξ, a, r)
  ν = stoichiometry(r)
  ξ .= ν * a
end

function propensity_derivatives!(b, x, r)
  for j in 1:size(b, 2)
    for k in eachindex(x)
      b[k, j] = compute_mass_action_deriv(x, r, j, k)
    end
  end
  return b
end

"""
  Gillespie 2001, Eq 26(a)
"""
function tauleap_DG2001(ξ, b, τ, a, ϵ)
  # τ .= ϵ .* a ./ abs.(b' * ξ)
  for j in eachindex(τ)
    τ[j] = zero(eltype(τ))
    for i in eachindex(ξ)
      τ[j] = τ[j] + ξ[i] * b[i, j]
    end
    τ[j] = ϵ * a / abs(τ[j])
  end
  return min(Inf, minimum(τ))
end
