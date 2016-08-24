"""
```
ODM(T)
```

Optimized Direct Method. Same as SSA, except the system is pre-simulated in order to sort reaction propensities from increasing to decreasing. This improves the search on the propensity vector when selecting the next reaction to fire.

### Arguments
- `T`: The simulation end time.

### Optional Arguments
- `init_steps`: Number of time steps to pre-simulate.
- `init_iters`: Number of iterations to simulate.
"""
type ODM <: ExactMethod
  # parameters
  end_time :: Float64
  nsteps   :: Int

  # state variables
  t     :: Float64
  steps :: Int

  # statistics
  avg_nsteps    :: Float64
  avg_stepsz :: Float64

  # metadata tags
  tags :: Vector{Symbol}

  function ODM(tf, nsteps)
    new(tf, nsteps, 0.0, 0, 0.0, 0.0, DEFAULT_EXACT)
  end
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
      # if any(x -> x < 0, Xt)
      #   error("t  = ",    get_time(algorithm),
      #         "\nXt = ", Xt,
      #         "\nμ  = ", μ,
      #         "\nV  = ", full(stoichiometry(r)),
      #         "\ndg = ", dependencies(r))
      # end
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

function sort!(r, ix)
  V  = stoichiometry(r)
  U  = coefficients(r)
  k  = scaled_rates(r)
  dg = dependencies(r)

  temp1 = V[:, ix]
  temp2 = U[:, ix]
  temp3 = dg[ix]

  for i in eachindex(V) # this is correct
    V[i] = temp1[i]
    U[i] = temp2[i]
  end

  for i in eachindex(dg) # this is not correct
    k[i], k[ix[i]], = k[ix[i]], k[i]
    dg[i] = temp3[i]
    for j in eachindex(dg[i])
      dg[i][j] = findfirst(ix, dg[i][j])
    end
  end

  return r
end
