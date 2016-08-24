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

function step!(algorithm::FRM, Xt::Vector, r::AbstractReactionSystem)
  a = propensities(r)

  if intensity(a) > 0
    τ, μ = select_reaction(algorithm, a)
    set_time!(algorithm, τ)

    if !done(algorithm)
      fire_reaction!(Xt, r, μ)
      update_dependent_propensities!(r, Xt, μ)
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
