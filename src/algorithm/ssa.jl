mutable struct SSA <: ExactMethod
  # parameters
  end_time :: Float64

  # state
  t      :: Float64

  # statistics
  stats_tracked :: Bool
  stats :: Dict{Symbol,Int}
  # metadata

  function SSA(end_time::AbstractFloat, stats_tracked)
    new(end_time, 0.0,
    stats_tracked,
      Dict{Symbol,Int}(
        :gillespie_steps => 0
    ))
  end
end

set_time!(algorithm::SSA, τ::AbstractFloat) = (algorithm.t = algorithm.t + τ)

##### implementation #####

function step!(algorithm::SSA, Xt::Vector, r::AbstractReactionSystem)
  a = propensities(r)

  if intensity(a) > 0
    τ = compute_stepsize(a)

    set_time!(algorithm, τ)

    if !isdone(algorithm)
      μ = select_reaction(a)
      fire_reaction!(Xt, r, μ)
      update_propensities!(r, Xt, μ)
    end

  elseif intensity(a) == 0
    algorithm.t = algorithm.end_time
  else
    throw(error("intensity = $(intensity(a)) < 0 at time $algorithm.t"))
  end

  if algorithm.stats_tracked
    algorithm.stats[:gillespie_steps] += 1
  end
  
  return nothing
end

##### next reaction #####

function compute_stepsize(a::PropensityVector)
  randexp() / intensity(a)
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
