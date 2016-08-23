"""
```
FRM(T)
```

Gillespie's First Reaction Method. Statistically equivalent to `SSA`, but slower computationally. The algorithm computes the time to the next reaction as the minimum the relative firing times for each reaction.

### Arguments
- `T`: The simulation end time.
"""
type FRM <: ExactMethod
  # parameters
  end_time :: Float64

  # state
  t      :: Float64
  nsteps :: Int

  # statistics
  avg_nsteps :: Float64
  avg_stepsz :: Float64

  # metadata
  tags :: Vector{Symbol}

  function FRM(T)
    new(T, 0.0, 0, 0.0, 0.0, DEFAULT_EXACT)
  end
end

set_time!(algorithm::FRM, τ::AbstractFloat) = (algorithm.t = algorithm.t + τ)

function init!(algorithm::FRM, a::PropensityVector)
  algorithm.t = 0.0
  return nothing
end

function step!(algorithm::FRM, Xt::Vector, r::AbstractReactionSystem)
  a = propensities(r)

  τ, μ = select_reaction(algorithm, a)
  set_time!(algorithm, τ)

  if !done(algorithm) && intensity(a) > 0
    fire_reaction!(Xt, r, μ)
    update_propensities!(r, Xt, μ)
  end

  return nothing
end

function select_reaction(::FRM, a::PropensityVector)
  min_val = Inf
  min_ind = 0

  for j in eachindex(a)
    temp = rand(Exponential(1 / a[j]))
    if temp < min_val
      min_val = temp
      min_ind = j
    end
  end

  return min_val, min_ind
end
