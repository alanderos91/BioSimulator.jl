type SSA <: ExactMethod
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

  function SSA(tf)
    new(tf, 0.0, 0, 0.0, 0.0, DEFAULT_EXACT)
  end
end

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
      update_propensities!(r, Xt, μ)
    end
  elseif intensity(a) == 0
    algorithm.t = algorithm.end_time
  else
    println("t = ", get_time(algorithm))
    println("a = ", a)
    println("Xt = ", Xt)
    error("intensity = ", intensity(a))
  end

  # update nsteps
  # update statistics

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

function chopdown(a::PropensityVector)
  jump = intensity(a) * rand()

  μ = length(a)

  while jump > 0
    jump -= a[μ]
    μ -= 1

    if μ == 0 && (jump > 0) println("a = ", a, " intensity = ", intensity(a)) end
  end

  return μ + 1
end
