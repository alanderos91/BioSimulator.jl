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
  nsteps       :: Int

  # statistics
  avg_nsteps    :: Float64
  avg_stepsz :: Float64

  # metadata tags
  tags :: Vector{Symbol}

  function NRM(end_time::AbstractFloat)
    new(end_time, 0.0, PriorityQueue(Int, Float64), 0, 0.0, 0.0, DEFAULT_EXACT)
  end
end

function NRM(;end_time=DEFAULT_TIME)
  return NRM(end_time)
end

get_reaction_times(algorithm::NRM) = algorithm.pq

set_time!(algorithm::NRM, t) = (algorithm.t = t)

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
  old_t = get_time(algorithm)
  if intensity(a) > 0
    pq = get_reaction_times(algorithm)

    μ, τ = peek(pq)

    # update algorithm variables
    set_time!(algorithm, τ)

    if !done(algorithm)
      fire_reaction!(Xt, r, μ)
      update_reaction_times!(algorithm, Xt, r, μ, τ)
    end

    # update nsteps
    nsteps!(algorithm)

    # update statistics
    compute_statistics!(algorithm, τ - old_t)

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
  t  = get_time(algorithm)

  dependents = dg[μ]

  for α in dependents
    temp = a[α]
    update_propensity!(r, Xt, α)

    if pq[α] < Inf
      pq[α] = t + (temp / a[α]) * (pq[α] - t)
    else
      pq[α] = τ + rand(Exponential(1 / a[α]))
    end
  end

  update_propensity!(r, Xt, μ)
  pq[μ] = τ + rand(Exponential(1 / a[μ]))

  return nothing
end
