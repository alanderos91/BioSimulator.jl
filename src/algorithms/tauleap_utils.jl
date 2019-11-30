abstract type TauLeapFormula end # <: Function

"""
Gillespie 2001, Eq 26(a)
"""
struct DG2001Eq26a{T1,T2} <: TauLeapFormula
  ξ::T1
  b::T2
  ϵ::Float64

  DG2001Eq26a(ξ::T1, b::T2, ϵ::Real) where {T1,T2} = new{T1,T2}(ξ, b, ϵ)
end

function DG2001Eq26a(number_species::Integer, number_jumps::Integer, ϵ::Real)
  ξ = zeros(number_species)
  b = zeros(number_species, number_jumps)

  return DG2001Eq26a(ξ, b, ϵ)
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
struct DGLP2003Eq6{U,B} <: TauLeapFormula
  V::U
  b::B
  ϵ::Float64
end

function DGLP2003Eq6(number_species, number_jumps, V, ϵ)
  b = zeros(number_species, number_jumps)

  return DGLP2003Eq6(V, b, ϵ)
end

function (formula::DGLP2003Eq6)(rates, total_rate)
  T = eltype(rates)
  V = formula.V
  b = formula.b
  ϵ = formula.ϵ
  τ = typemax(T)

  A1 = ϵ * total_rate
  A2 = A1^2

  for j in eachindex(rates)
    μ, σ² = zero(T), zero(T)
    f = f_calculation(V, b, j)
    for k in eachindex(rates)
      μ  = μ + f * rates[k]
      σ² = σ² + f^2 * rates[k]
    end
    τ′ = min(A1 / abs(μ), A2 / σ²)
    τ  = min(τ, τ′)
  end

  return τ
end

function f_calculation(V::DenseMatrix{Int}, b, j)
  f = zero(eltype(b))

  @simd for i in 1:size(b, 1)
    @inbounds @fastmath f = f + b[i, j] * V[i, j]
  end

  return f
end

function f_calculation(V::SparseMatrixCSC{Int,Int}, b, j)
  f = zero(eltype(b))
  rv = rowvals(V)
  nz = nonzeros(V)

  for i in nzrange(V, j)
    @inbounds @fastmath f = f + b[rv[i], j] * nz[i]
  end

  return f
end

function update!(formula::DGLP2003Eq6, state, model, stoichiometry, rates, total_rate)
  b = formula.b

  n, m = size(b)

  for j in 1:m
    for i in 1:n
      b[i,j] = rate_derivative(model, state, i, j)
    end
  end

  return nothing
end

"""
Unknown formula?
"""
struct GenericLeapFormula
  k::Vector{Float64}
  U::SparseMatrixCSC{Int,Int}
  V::SparseMatrixCSC{Int,Int}
  dxdt::Vector{Float64}
  drdt::Vector{Float64}
  ϵ::Float64
end

# function GenericLeapFormula(number_species::Integer, number_jumps::Integer, U, V, ϵ)
#   dxdt = zeros(number_species)
#   drdt = zeros(number_jumps)

#   return GenericLeapFormula(U, V, dxdt, drdt, ϵ)
# end

function (formula::GenericLeapFormula)(rates, total_rate)
  drdt = formula.drdt
  k = formula.k
  τ = typemax(total_rate)
  ϵ = formula.ϵ

  for j in eachindex(rates)
    A = ϵ * max(rates[j], k[j])
    B = abs(drdt[j])

    τ = min(τ, A / B)
  end

  return τ
end

function update!(formula::GenericLeapFormula, state, model, V, rates, total_rate)
  dxdt = formula.dxdt
  drdt = formula.drdt
  # V = formula.V
  U = formula.U

  update_mean_derivatives!(dxdt, V, rates)
  update_time_derivatives!(drdt, state, model, dxdt)
end

update_mean_derivatives!(dxdt, V, rates) = (dxdt .= V * rates)

function update_time_derivatives!(drdt, state, model, dxdt)
  for j in eachindex(drdt)
    drdt[j]  = 0.0
    for i in eachindex(state)
      drdx = rate_derivative(model, state, i, j)
      drdt[j] = drdt[j] + drdx * dxdt[i]
    end
  end

  return drdt
end

##### post-leap check for valid step #####

struct RejectionThinning{T1,T2,F1,F2,F3,F4}
  threshold::T1
  proposal::T2
  is_accepted::F1
  is_invalid::F2
  execute_leap!::F3
  reverse_leap!::F4
end

function RejectionThinning(threshold::Real, proposal, execute_leap!, reverse_leap!)
  is_accepted = Base.Fix2(<, threshold)     # x -> x < threshold
  is_negative = Base.Fix2(<, 0)             # x -> x < 0
  is_invalid  = Base.Fix1(any, is_negative) # x -> any(is_negative, x)

  RejectionThinning(threshold, proposal, is_accepted, is_invalid, execute_leap!, reverse_leap!)
end

function (f::RejectionThinning)(v::AbstractVector, s::Real)
  # unpack
  threshold = f.threshold
  is_invalid = f.is_invalid
  is_accepted = f.is_accepted
  execute_leap! = f.execute_leap!
  reverse_leap! = f.reverse_leap!
  proposal = f.proposal

  execute_leap!(proposal, v)

  while is_invalid(proposal)
    reverse_leap!(proposal, v)

    # reject events in the proposed leap
    for j in eachindex(v)
      # generate a stream of iid uniform deviates for the rejection step
      event = (rand() for j in 1:v[j])

      # count the number of accepted events
      v[j] = count(is_accepted, event)
    end

    # contract the leap by the threshold parameter
    s = s * threshold

    # apply new update
    execute_leap!(proposal, v)
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

function extract_coefficients(model::ReactionSystem)
  reaction = model.reactions

  number_entries = sum(length(reaction[j].reactants) for j in eachindex(reaction))

  I = zeros(Int, number_entries) # species index
  J = zeros(Int, number_entries) # reaction index
  V = zeros(Int, number_entries) # net stoichiometry

  idx = 1

  for j in eachindex(reaction)
    for (i, v) in reaction[j].reactants

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

function forward_leap!(x, V, n)
  x .+= V*n
end

function backward_leap!(x, V, n)
  x .-= V*n
end
