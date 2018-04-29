"""
```
FRM
```

Gillespie's First Reaction Method.

Statistically equivalent to `SSA`, but more computationally expensive: it computes the time to the next reaction as the minimum waiting time relative to the next firing times of each reaction.

### Internals
- `end_time`: The termination time, supplied by a user.
- `t`: The current simulation time.
"""
mutable struct FRM <: ExactMethod
  # parameters
  end_time :: Float64

  # state
  t :: Float64

  # statistics
  stats :: Dict{Symbol,Int}

  function FRM(end_time::AbstractFloat)
    new(end_time, 0.0,
      Dict{Symbol,Int}(
        :gillespie_steps => 0
    ))
  end
end

FRM(end_time; na...) = FRM(end_time)

set_time!(algorithm::FRM, τ::AbstractFloat) = (algorithm.t = algorithm.t + τ)

function step!(algorithm::FRM, Xt::Vector, r::AbstractReactionSystem)
  a = propensities(r)

  if intensity(a) > 0
    τ, μ = select_reaction(algorithm, a)
    set_time!(algorithm, τ)

    if !done(algorithm)
      fire_reaction!(Xt, r, μ)
      update_propensities!(a, r, Xt, μ)
    end

  elseif intensity(a) == 0
    algorithm.t = algorithm.end_time
  else
    throw(Error("intensity = $(intensity(a)) < 0 at time $algorithm.t"))
  end

  algorithm.stats[:gillespie_steps] += 1
  
  return nothing
end

@inbounds function select_reaction(::FRM, a::PropensityVector)
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
