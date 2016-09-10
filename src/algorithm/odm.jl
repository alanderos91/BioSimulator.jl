"""
```
ODM(end_time=1.0, n_steps=1000)
```

Optimized Direct Method. Same as `SSA`, except the system is presimulated to sort reaction propensities from increasing to decreasing. This improves the search on the CMF when selecting the next reaction to fire.

### Arguments
- `end_time`: The simulation end time.
- `n_steps`: Number of time steps to presimulate.
"""
type ODM <: ExactMethod
  # parameters
  end_time :: Float64
  n_steps   :: Int

  # state variables
  t      :: Float64
  nsteps :: Int

  # statistics
  avg_nsteps :: Mean{EqualWeight}

  # metadata tags
  tags :: Vector{Symbol}

  function ODM(end_time::AbstractFloat, n_steps::Integer)
    new(end_time, n_steps, 0.0, 0, Mean(), DEFAULT_EXACT)
  end
end

function ODM(;end_time=DEFAULT_TIME, n_steps=1000)
  return ODM(end_time, n_steps)
end

set_time!(algorithm::ODM, τ) = (algorithm.t = algorithm.t + τ)

function init!(algorithm::ODM, Xt, r)
  nsteps = algorithm.nsteps
  reaction_events = zeros(Int, size(stoichiometry(r), 2))

  presimulate!(reaction_events, Xt, r, nsteps)

  ix = sortperm(reaction_events)
  sort!(r, ix)

  return nothing
end

function step!(algorithm::ODM, Xt, r)
  a = propensities(r)

  if intensity(a) > 0

    τ  = compute_stepsize(a)

    # update algorithm variables
    set_time!(algorithm, τ)

    if !done(algorithm)
      μ = select_reaction(a)
      fire_reaction!(Xt, r, μ)
      update_dependent_propensities!(r, Xt, μ)
    end

    algorithm.nsteps += 1

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

function presimulate!(
  reaction_events :: Vector{Int},
  Xt              :: Vector{Int},
  r               :: AbstractReactionSystem,
  nsteps          :: Integer
  )

  update_all_propensities!(r, Xt)

  a = propensities(r)

  for i = 1:nsteps
    if intensity(a) > 0
      μ = select_reaction(a)
      fire_reaction!(Xt, r, μ)
      update_dependent_propensities!(r, Xt, μ)

      if islossy(a)
        a.intensity = sum(a)
        a.error_bound = zero(eltype(a))
      end

      reaction_events[μ] = reaction_events[μ] + 1
    else
      break
    end
  end

  return reaction_events
end

@inbounds function sort!(r, ix)
  V  = stoichiometry(r)
  U  = coefficients(r)
  k  = scaled_rates(r)
  dg = dependencies(r)

  temp1 = V[:, ix]
  temp2 = U[:, ix]
  temp3 = k[ix]
  temp4 = dg[ix]

  for i in eachindex(V) # this is correct
    V[i] = temp1[i]
    U[i] = temp2[i]
  end

  for i in eachindex(dg) # this is not correct
    k[i]  = temp3[i]
    dg[i] = temp4[i]
    for j in eachindex(dg[i])
      dg[i][j] = findfirst(ix, dg[i][j])
    end
  end

  return r
end
