"""
```
SSA
```

Gillespie's Direct Method (SSA). Simulates a system of coupled reactions and species by computing the time to the next reaction and searching on the CMF.

### Internals
- `end_time`: The termination time, supplied by a user.
- `t`: The current simulation time.
"""
type SSA <: ExactMethod
  # parameters
  end_time :: Float64

  # state
  t      :: Float64

  # statistics

  # metadata

  function SSA(end_time::AbstractFloat)
    new(end_time, 0.0)
  end
end

SSA(end_time; na...) = SSA(end_time)

set_time!(algorithm::SSA, τ::AbstractFloat) = (algorithm.t = algorithm.t + τ)

##### implementation #####

function step!(algorithm::SSA, Xt::Vector, r::AbstractReactionSystem)
  a = propensities(r)

  if intensity(a) > 0
    τ = compute_stepsize(a)

    set_time!(algorithm, τ)

    if !done(algorithm)
      μ = select_reaction(a)
      fire_reaction!(Xt, r, μ)
      update_propensities!(a, r, Xt, μ)
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

##### next reaction #####

function compute_stepsize(a::PropensityVector)
  rand(Exponential(1 / intensity(a)))
end

##### selecting reaction #####

function select_reaction(a::PropensityVector)
  chopdown(a)
end

function buildup(a::PropensityVector)
  jump = intensity(a) * rand()
  asum = zero(eltype(a))

  μ = 1

  while asum < jump
    asum += a[μ]
    μ += 1
  end

  return μ - 1
end

@inbounds @fastmath function chopdown(a::PropensityVector)
  jump = intensity(a) * rand()

  μ = length(a)

  while jump > 0
    jump -= a[μ]
    μ -= 1
  end

  return μ + 1
end
