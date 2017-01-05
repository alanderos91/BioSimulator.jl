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

  # state variables
  t      :: Float64

  # statistics

  function ODM(end_time::AbstractFloat)
    new(end_time, 0.0)
  end
end

function ODM(;end_time=0.0, na...)
  if end_time == 0.0
    error("end_time argument must be positive.")
  end
  return ODM(end_time)
end

set_time!(algorithm::ODM, τ) = (algorithm.t = algorithm.t + τ)

function init!(algorithm::ODM, Xt, r)
  reaction_events = zeros(Int, size(stoichiometry(r), 2))

  presimulate!(reaction_events, Xt, r, end_time(algorithm))

  ix = sortperm(reaction_events)
  _sort!(r, ix)

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

function presimulate!(
  reaction_events :: Vector{Int},
  Xt              :: Vector{Int},
  r               :: AbstractReactionSystem,
  end_time        :: AbstractFloat
  )

  a = propensities(r)
  t = zero(typeof(end_time))

  update_all_propensities!(a, r, Xt)

  while t < end_time
    if intensity(a) > 0
      τ = compute_stepsize(a)
      t = t + τ

      if t < end_time
        μ = select_reaction(a)
        fire_reaction!(Xt, r, μ)
        update_propensities!(a, r, Xt, μ)
        reaction_events[μ] = reaction_events[μ] + 1
      end
    else
      break
    end
  end

  return reaction_events
end

@inbounds function _sort!(r, ix)
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
