mutable struct FRM <: ExactMethod
  # parameters
  end_time :: Float64

  # state
  t :: Float64

  # statistics
  stats_tracked :: Bool
  stats :: Dict{Symbol,Int}

  function FRM(end_time::AbstractFloat, stats_tracked)
    new(end_time, 0.0,
      stats_tracked,
      Dict{Symbol,Int}(
        :gillespie_steps => 0
    ))
  end
end

set_time!(algorithm::FRM, τ::AbstractFloat) = (algorithm.t = algorithm.t + τ)

function step!(algorithm::FRM, Xt::Vector, r::AbstractReactionSystem)
  a = propensities(r)

  if intensity(a) > 0
    τ, μ = select_reaction(algorithm, a)
    set_time!(algorithm, τ)

    if !isdone(algorithm)
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

@inbounds function select_reaction(::FRM, a::PropensityVector)
  min_val = Inf
  min_ind = 0

  for j in eachindex(a)
    temp = randexp() / a[j]
    if temp < min_val
      min_val = temp
      min_ind = j
    end
  end

  return min_val, min_ind
end
