abstract type TauLeapFormula end

"""
Gillespie 2001, Eq 26(a)
"""
struct DG2001Eq26a{T1,T2} <: TauLeapFormula
  ξ::T1
  b::T2
  ϵ::Float64
end

function DG2001Eq26a(number_species, number_jumps)
  ξ = zeros(number_species)
  b = zeros(number_species, number_jumps)

  return DG2001Eq26a(ξ, b)
end

function (formula::DG2001Eq26a)(rates, total_rate)
  ξ = formula.ξ
  b = formula.b
  ϵ = formula.ϵ

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

function update!(formula::DG2001Eq26a, state, model, stoichiometry, rates, total_rate)
  ξ = formula.ξ
  b = formula.b

  mul!(ξ, stoichiometry, rates)

  n, m = size(b)

  for j in 1:m
    for i in 1:n
      b[i,j] = rate_derivative(model, state, i, j)
    end
  end

  return nothing
end

"""
  Gillespie and Petzold 2003, Eq 6
"""
struct DGLP2003{V,B} <: TauLeapCache
  ν::V
  b::B
end

function tauleap(cache::DGLP2003, rates, total_rate, ϵ)
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
      μ  = μ + f * rates[k]
      σ² = σ² + f^2 * rates[k]
    end
    τ′ = min(A1 / abs(μ), A2 / σ²)
    τ  = min(τ, τ′)
  end

  return τ
end

function f_calculation(ν::DenseMatrix{Int}, b, j)
  f = zero(eltype(b))

  @simd for i in 1:size(b, 1)
    @inbounds @fastmath f = f + b[i, j] * ν[i, j]
  end

  return f
end

function f_calculation(ν::SparseMatrixCSC{Int,Int}, b, j)
  f = zero(eltype(b))
  rv = rowvals(ν)
  nz = nonzeros(ν)

  for i in nzrange(ν, j)
    @inbounds @fastmath f = f + b[rv[i], j] * nz[i]
  end

  return f
end

##### post-leap check for valid step #####

# TODO
function is_bad_tauleap(proposal)
  return any(Base.Fix2(<,0), proposal)
end

struct RejectionThinning{T,F1,F2}
  threshold::T
  is_accepted::F1
  is_invalid::F2
end

function RejectionThinning(threshold::Real)
  is_accepted = Base.Fix2(<, threshold)     # x -> x < threshold
  is_negative = Base.Fix2(<, 0)             # x -> x < 0
  is_invalid  = Base.Fix1(any, is_negative) # x -> any(is_negative, x)

  RejectionThinning(threshold, is_accepted, is_invalid)
end

function (f::RejectionThinning)(x, v::AbstractVector, s::Real)
  is_invalid = f.is_invalid
  is_accepted = f.is_accepted

  while is_invalid(x)
    # reject events in the proposed leap
    for j in eachindex(v)
      itr = (rand() for j in 1:v[j])
      v[j] = count(is_accepted, itr)
    end

    # contract the leap by the threshold parameter
    s = s * f.threshold

    # propose update
    mul!(x, A, v)
  end

  return v, s
end

##### extracting stoichiometry matrix #####

function extract_net_stoichiometry(model::ReactionSystem)
  reaction = model.reactions

  number_entries = sum(length(reaction[j].net_change) for j in eachindex(reaction))

  I = zeros(Int, number_entries) # species index
  J = zeros(Int, number_entries) # reaction index
  V = zeros(Int, number_entries) # net stoichiometry

  idx = 1

  for j in eachindex(reaction)
    for (i, v) in reaction[j].net_change

      I[idx] = i
      J[idx] = j
      V[idx] = v

      idx += 1
    end
  end

  return sparse(I, J, V)
end

##### applying leap updates #####

"""
`ApplyLeapUpdate`

A type-stable closure used to update system state.
Specifically, this type is used to cache the net stoichiometry matrix.

##### Usage

```
execute_leap! = ApplyLeapUpdate(mul!, V)
execute_leap!(x, n) # mul!(x, V, n)
```
"""
struct ApplyLeapUpdate{F,T} <: Function
  f::F
  V::T

  ApplyLeapUpdate(f::F, V::T) where {F,T} = new{F,T}(f, V)
  ApplyLeapUpdate(f::Type{F}, x::T) where {F,T} = new{Type{F},T}(f, V)
end

(g::ApplyLeapUpdate)(x, n) = g.f(x, g.V, n)
