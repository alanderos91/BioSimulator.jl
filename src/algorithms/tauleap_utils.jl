abstract type TauLeapCache end

"""
Gillespie 2001, Eq 26(a)
"""
struct DG2001{X,B} <: TauLeapCache
  ξ::X
  b::B
end

function tauleap(cache::DG2001, τ, rates, ϵ)
  ξ = cache.ξ
  b = cache.b
  τ = typemax(eltype(rates))

  for j in eachindex(rates)
    proposal = zero(τ)
    for i in eachindex(ξ)
      proposal = proposal + ξ[i] * b[i, j]
    end
    proposal = ϵ * total_rate / abs(proposal)

    if proposal < τ
      τ = proposal
    end
  end

  return τ
end

"""
  Gillespie and Petzold 2003, Eq 6
"""
struct DGLP2003{V,B} <: TauLeapCache
  ν::V
  b::B
end

function tauleap(cache::DGLP2003)
  T = eltype(rates)
  A1 = ϵ * total_rate
  A2 = A1^2
  
  ν = cache.ν
  b = cache.b
  τ = typemax(T)

  for j in eachindex(rates)
    μ, σ² = zero(T), zero(T)
    f = f_calculation(ν, b, j)
    for k in eachindex(rates)
      μ  = μ + f * a[k]
      σ² = σ² + f^2 * a[k]
    end
    τ′ = min(A1 / abs(μ), A2 / σ²)
    τ  = min(τ, τ′)
  end
  return τ
end