import Base.Collections: PriorityQueue, peek
import Base.Order: ForwardOrdering

"""
```
NRM(end_time=1.0)
```

Gibson and Bruck's Next Reaction Method, statistically equivalent to the original `SSA`. Provides better computational efficiency on networks with loosely connected reactions.

### Arguments
- `end_time`: The end time for the simulation.
"""
type NRM <: ExactMethod
  # parameters
  end_time :: Float64

  # state variables
  t            :: Float64
  pq           :: PriorityQueue{Int,Float64,ForwardOrdering}

  # statistics

  function NRM(end_time::AbstractFloat)
    new(end_time, 0.0, PriorityQueue(Int, Float64))
  end
end

function NRM(;end_time=0.0, na...)
  if end_time == 0.0
    error("end_time argument must be positive.")
  end
  return NRM(end_time)
end

get_reaction_times(algorithm::NRM) = algorithm.pq

set_time!(algorithm::NRM, t) = (algorithm.t = t)

function init!(algorithm::NRM, Xt, r)
  dg = dependencies(r)
  a  = propensities(r)
  pq = get_reaction_times(algorithm)

  for j in eachindex(a)
    pq[j] = zero(eltype(a))
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
    println("t = ", get_time(algorithm))
    println("a = ", a)
    println("Xt = ", Xt)
    error("intensity = ", intensity(a))
  end

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

  return nothing
end
