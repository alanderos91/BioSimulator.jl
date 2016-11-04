"""
```
simulate{T<:Algorithm}(model::Network, algorithm::Type{T};
           time::AbstractFloat=1.0,
           epochs::Integer=1,
           trials::Integer=1,
           algvars...)
```
### Arguments
- `model`: The `Network` to simulate.
- `algorithm`: The algorithm used to carry out simulation. One of `SSA`, `FRM`, `NRM`, `ODM`, or `SAL`.

### Optional Arguments
- `time`: The amount of time to simulate the model, in units of the model.
- `epochs`: The number of times to sample the vector of counts.
- `trials`: The number of independent realizations to generate.
- `kwargs`: Additional keyword arguments specific to each algorithm.

"""
function simulate{T}(model::Network, algorithm::Type{T}=SSA;
           time::AbstractFloat=1.0,
           epochs::Integer=1,
           trials::Integer=1,
           kwargs...)

    # extract model information
    c = n_species(model)
    d = n_reactions(model)

    species = species_list(model)
    reactions = reaction_list(model)

    # create simulation data structures
    x0, rxn, output, id2ind = make_datastructs(species, reactions, c, d)

    # initialize algorithm
    alg = algorithm(;end_time=time, kwargs...)

    # run simulation
    if nworkers() == 1
      output = PartialHistory(DenseArray, c, epochs+1, trials, 0.0, time, id2ind)
      serial_simulate(output, x0, alg, deepcopy(x0), rxn)
    else
      output = PartialHistory(SharedArray, c, epochs+1, trials, 0.0, time, id2ind)
      parallel_simulate(output, x0, alg, deepcopy(x0), rxn)
    end
end

function make_datastructs(species, reactions, c, d)
  # state vector
  x0, id, id2ind = make_species_vector(species)

  # reactions
  if d <= 8
    rxn = DenseReactionSystem(reactions, id2ind, c, d)
  else
    rxn = SparseReactionSystem(reactions, id2ind, c, d)
  end

  return x0, rxn, id, id2ind
end

function serial_simulate(
  output    :: PartialHistory,
  Xt        :: Vector{Int},
  algorithm :: Algorithm,
  X0        :: Vector{Int},
  r         :: AbstractReactionSystem) # nrlz is encoded in PartialHistory; refactor

  a = propensities(r)

  init!(algorithm, Xt, r)
  for i in 1:size(output.data, 3)
    # setup
    copy!(Xt, X0)
    update_all_propensities!(a, r, Xt)
    reset!(algorithm, a)
    interval = 1

    while !done(algorithm)
      interval = update!(output, Xt, get_time(algorithm), interval, i)
      step!(algorithm, Xt, r)
    end

    interval = update!(output, Xt, get_time(algorithm), interval, i)
  end

  return output
end

function parallel_simulate(
  output    :: PartialHistory,
  Xt        :: Vector{Int},
  algorithm :: Algorithm,
  X0        :: Vector{Int},
  r         :: AbstractReactionSystem) # nrlz is encoded in PartialHistory; refactor

  init!(algorithm, Xt, r)

  @sync for pid in procs(output.data)
      @async remotecall_fetch(pid, trajectory_shared_chunk!, output, Xt, algorithm, X0, r)
  end

  return output
end

# Adapted from SharedArrays documentation
function nrlz_partition(output::PartialHistory)
    q = output.data
    idx = indexpids(q)
    if idx == 0
        # This worker is not assigned a piece
        return 1:0, 1:0
    end
    nchunks = length(procs(q))
    splits = [round(Int, s) for s in linspace(0,size(q,3),nchunks+1)]
    splits[idx]+1:splits[idx+1]
end

function trajectory_chunk!(
  output :: PartialHistory,
  Xt        :: Vector{Int},
  algorithm :: Algorithm,
  X0        :: Vector{Int},
  r         :: AbstractReactionSystem,
  krange    :: UnitRange
  )

  a = propensities(r)

  for i in krange
    copy!(Xt, X0)
    update_all_propensities!(a, r, Xt)
    reset!(algorithm, a)
    interval = 1

    while !done(algorithm)
      interval = update!(output, Xt, get_time(algorithm), interval, i)
      step!(algorithm, Xt, r)
    end

    interval = update!(output, Xt, get_time(algorithm), interval, i)
  end

  return output
end

# the wrapper
@inline function trajectory_shared_chunk!(
  output :: PartialHistory,
  Xt        :: Vector{Int},
  algorithm :: Algorithm,
  X0        :: Vector{Int},
  r         :: AbstractReactionSystem
  )
  trajectory_chunk!(output, Xt, algorithm, X0, r, nrlz_partition(output))
end
