mutable struct OTL <: TauLeapMethod
  # parameters
  end_time :: Float64
  ϵ :: Float64
  δ :: Float64
  β :: Float64

  # state variables
  t :: Float64

  # cache
  # ξ :: Vector{Float64}
  b :: Matrix{Float64}
  # τ :: Vector{Float64}
  events :: Vector{Int}

  # statistics
  stats_tracked :: Bool
  stats :: Dict{Symbol,Int}

  # function OTL(end_time, ϵ, δ, β, stats_tracked)
  #   new(end_time, ϵ, δ, β,
  #   0.0, Float64[], Matrix{Float64}(0,0), Float64[], Int[],
  #   stats_tracked,
  #   Dict{Symbol,Int}(
  #     :negative_excursions => 0,
  #     :contractions => 0,
  #     :leaping_steps => 0,
  #     :gillespie_steps => 0
  #   ))
  # end
  function OTL(end_time, ϵ, δ, β, stats_tracked)
    new(end_time, ϵ, δ, β,
    0.0, Matrix{Float64}(undef,0,0), Int[],
    stats_tracked,
    Dict{Symbol,Int}(
      :negative_excursions => 0,
      :contractions => 0,
      :leaping_steps => 0,
      :gillespie_steps => 0
    ))
  end
end

set_time!(algorithm::OTL, τ) = (algorithm.t = algorithm.t + τ)

function init!(x :: OTL, Xt, r)
  c = length(Xt)
  d = size(stoichiometry(r), 2)

  # x.ξ = zeros(c)    # store expected change
  x.b = zeros(c, d) # store partial derivatives
  # x.τ = zeros(d)    # compute leap size
  x.events = zeros(Int, d)

  return nothing
end

function reset!(algorithm::OTL, Xt, r)
  update_all_propensities!(r, Xt)
  propensity_derivatives!(algorithm.b, Xt, r)
  a = propensities(r)
  algorithm.t = zero(algorithm.t)

  return nothing
end

function step!(algorithm :: OTL, Xt, r)
  a = propensities(r)
  iscritical = (intensity(a) < algorithm.δ)

  # compute leap size and check for end of interval
  τ = tauleap_DGLP2003(stoichiometry(r), algorithm.b, a.cache, a.intensity, algorithm.ϵ)
  τ = min(τ, end_time(algorithm) - get_time(algorithm))

  if intensity(a) > 0
    if τ < algorithm.δ / intensity(a)
      # Gillespie update
      if algorithm.stats_tracked
        algorithm.stats[:gillespie_steps] += 1
      end
      τ = randexp() / intensity(a)
      set_time!(algorithm, τ)

      if !done(algorithm)
        μ = select_reaction(a)
        fire_reaction!(Xt, r, μ)
        update_propensities!(r, Xt, μ)
      end
    else
      # τ-leap update
      if algorithm.stats_tracked
        algorithm.stats[:leaping_steps] += 1
      end

      # compute jumps
      for j in eachindex(algorithm.events)
        algorithm.events[j] = poisrand(τ * a[j])
      end

      isbadleap = is_badleap(Xt, r, algorithm.events)
      if isbadleap && algorithm.stats_tracked
        algorithm.stats[:negative_excursions] += 1
      end
    
      while isbadleap
        if algorithm.stats_tracked
          algorithm.stats[:contractions] += 1
        end
        contract!(algorithm.events, algorithm.β)
        τ = τ * algorithm.β
        isbadleap = is_badleap(Xt, r, algorithm.events)
      end
    
      fire_reactions!(Xt, r, algorithm.events)

      # state dependent updates
      update_all_propensities!(r, Xt)
      propensity_derivatives!(algorithm.b, Xt, r)

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
      @inbounds b[k, j] = compute_mass_action_deriv(x, r, j, k)
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

"""
  Gillespie and Petzold 2003, Eq 6
"""
function tauleap_DGLP2003(ν, b, a, a0, ϵ)
  T = eltype(a0)
  
  A1 = ϵ * a0
  A2 = A1^2
  τ = typemax(T)

  for j in eachindex(a)
    μ, σ² = zero(T), zero(T)
    for k in eachindex(a)
      f = f_calculation(ν, b, j)

      @inbounds @fastmath μ  = μ + f * a[k]
      @inbounds @fastmath σ² = σ² + f^2 * a[k]
    end
    τ′ = min(A1 / abs(μ), A2 / σ²)
    τ  = min(τ, τ′)
  end
  return τ
end

function f_calculation(ν :: DenseMatrix{Int}, b, j)
  f = zero(eltype(b))

  @simd for i in 1:size(b, 1)
    @inbounds @fastmath f = f + b[i, j] * ν[i, j]
  end

  return f
end

function f_calculation(ν :: SparseMatrixCSC{Int,Int}, b, j)
  f = zero(eltype(b))
  rv = rowvals(ν)
  nz = nonzeros(ν)

  for i in nzrange(ν, j)
    @inbounds @fastmath f = f + b[rv[i], j] * nz[i]
  end

  return f
end
