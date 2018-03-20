"""
```
NRM
```

Gibson and Bruck's Next Reaction Method, statistically equivalent to `SSA`. It provides better computational efficiency on networks with loosely connected reactions.

### Internals
- `end_time`: The termination time, supplied by a user.
- `t`: The current simulation time.
- `pq`: A priority queue that sorts reaction according to the their next firing times.
"""
mutable struct NRM <: ExactMethod
  # parameters
  end_time :: Float64

  # state variables
  t        :: Float64
  pq       :: PriorityQueue{Int,Float64,Base.Order.ForwardOrdering}

  # statistics
  stats :: Dict{Symbol,Int}

  function NRM(end_time::AbstractFloat)
    new(end_time, 0.0, PriorityQueue{Int, Float64}(),
      Dict{Symbol,Int}(
        :gillespie_steps => 0
    ))
  end
end

NRM(end_time; na...) = NRM(end_time)

get_reaction_times(algorithm::NRM) = algorithm.pq

set_time!(algorithm::NRM, t) = (algorithm.t = t)

function init!(algorithm::NRM, Xt, r)
  dg = dependencies(r)
  a  = propensities(r)
  pq = get_reaction_times(algorithm)

  for j in eachindex(a)
    pq[j] = zero(eltype(a))
    # This is a check for reactions with constant propensities. We do not need to update their propensities, but we do need to update their putative times.
    if j ∉ dg[j]
      push!(dg[j], j)
    end
  end
end

function reset!(algorithm::NRM, a::PVec)
  algorithm.t = 0.0
  pq = algorithm.pq

  for j in eachindex(a)
    pq[j] = rand(Exponential(1 / a[j]))
  end

  return nothing
end

function step!(algorithm::NRM, Xt::Vector, r::AbstractReactionSystem)
  a = propensities(r)

  if intensity(a) > 0
    pq = get_reaction_times(algorithm)

    μ, τ = peek(pq)

    # update algorithm variables
    set_time!(algorithm, τ)

    if !done(algorithm)
      fire_reaction!(Xt, r, μ)
      update_reaction_times!(algorithm, Xt, r, μ, τ)
    end

  elseif intensity(a) == 0
    algorithm.t = algorithm.end_time
  else
    throw(Error("intensity = $(intensity(a)) < 0 at time $algorithm.t"))
  end

  algorithm.stats[:gillespie_steps] += 1
  
  return nothing
end

function update_reaction_times!(algorithm::NRM, Xt, r, μ, τ)
  a  = propensities(r)
  k  = scaled_rates(r)
  dg = dependencies(r)
  pq = get_reaction_times(algorithm)

  dependents = dg[μ]
  T = eltype(a)

  for α in dependents
    oldval = a[α]
    old_t  = pq[α]
    update_propensity!(a, r, Xt, α)

    if α != μ && oldval != zero(T)

      if a[α] > zero(T)
        pq[α] = τ + (oldval / a[α]) * (old_t - τ)
      else
        pq[α] = Inf
      end

    else

      if a[α] > zero(T)
        pq[α] = τ + rand(Exponential(1 / a[α]))
      else
        pq[α] = Inf
      end

    end
  end

  !isstable(a) && update_all_propensities!(a, r, Xt)

  return nothing
end
