import Base.Collections: PriorityQueue, peek
import Base.Order: ForwardOrdering

"""
```
NRM(T)
```

Gibson and Bruck's Next Reaction Method, equivalent to the original `SSA`.

### Arguments
- `T`: The simulation end time.
"""
type NRM <: ExactMethod
  # parameters
  end_time :: Float64

  # state variables
  t            :: Float64
  pq           :: PriorityQueue{Int,Float64,ForwardOrdering}
  steps        :: Int

  # statistics
  avg_nsteps    :: Float64
  avg_stepsz :: Float64

  # metadata tags
  tags :: Vector{Symbol}

  function NRM(tf)
    new(tf, 0.0, PriorityQueue(Int, Float64), 0, 0.0, 0.0, DEFAULT_EXACT)
  end
end

get_reaction_times(algorithm::NRM) = algorithm.pq

set_time!(algorithm::NRM, t) = (algorithm.t = t)

function init!(algorithm::NRM, a::PropensityVector)
  algorithm.t = 0.0
  pq = algorithm.pq

  for j in eachindex(a)
    pq[j] = rand(Exponential(1 / a[j]))
  end

  return nothing
end

function step!(algorithm::NRM, Xt::Vector, r::AbstractReactionSystem)
  a = propensities(r)

  pq = get_reaction_times(algorithm)

  μ, τ = peek(pq)

  # update algorithm variables
  set_time!(algorithm, τ)

  if !done(algorithm) && intensity(a) > 0
    fire_reaction!(Xt, r, μ)
    update_reaction_times!(algorithm, Xt, r, μ, τ)
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
    a[α] = compute_mass_action(Xt, r, α)

    if pq[α] < Inf
      pq[α] = t + (temp / a[α]) * (pq[α] - t)
    else
      pq[α] = τ + rand(Exponential(1 / a[α]))
    end
  end

  a[μ]  = compute_mass_action(Xt, r, μ)
  pq[μ] = τ + rand(Exponential(1 / a[μ]))

  return nothing
end
