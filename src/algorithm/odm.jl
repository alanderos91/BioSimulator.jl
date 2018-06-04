"""
```
ODM
```

Optimized Direct Method. Similar to `SSA`, with the added benefit of sorting reactions according to their propensities over time. This improves the search on the CMF when selecting the next reaction to fire.

### Internals
- `end_time`: The termination time, supplied by a user.
- `t`: The current simulation time.
"""
mutable struct ODM <: ExactMethod
  # parameters
  end_time :: Float64

  # state variables
  t      :: Float64

  # statistics
  stats_tracked :: Bool
  stats :: Dict{Symbol,Int}

  function ODM(end_time::AbstractFloat, stats_tracked)
    new(end_time, 0.0,
        stats_tracked,
        Dict{Symbol,Int}(
        :gillespie_steps => 0
    ))
  end
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
    throw(Error("intensity = $(intensity(a)) < 0 at time $algorithm.t"))
  end

  if algorithm.stats_tracked
    algorithm.stats[:gillespie_steps] += 1
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
  V = stoichiometry(r)[:, ix]
  U = coefficients(r)[:, ix]
  k = scaled_rates(r)[ix]
  g = dependencies(r)[ix]

  for i in eachindex(V)
    r.stoichiometry[i] = V[i]
    r.coefficients[i] = U[i]
  end

  for i in eachindex(ix)
    r.scaled_rates[i] = k[i]
    r.dependencies[i] = g[i]
    for j in eachindex(g[i])
      r.dependencies[i][j] = findfirst(ix, g[i][j])
    end
  end

  if issparse(V)
    dropzeros!(r.stoichiometry)
    dropzeros!(r.coefficients)
  end

  return r
end
