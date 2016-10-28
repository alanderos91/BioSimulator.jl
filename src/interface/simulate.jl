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
- `algvars`: Additional keyword arguments specific to each algorithm.

"""
function simulate{T}(model::Network, algorithm::Type{T}=SSA;
           time::AbstractFloat=1.0,
           epochs::Integer=1,
           trials::Integer=1,
           algvars...)

    # extract model information
    c = n_species(model)
    d = n_reactions(model)

    species = species_list(model)
    reactions = reaction_list(model)

    # create simulation data structures
    x0, rxn, output = make_datastructs(species, reactions, c, d)

    # initialize algorithm
    alg = call(algorithm, algvars...)

    # run simulation
    simulate(...)
end

function make_datastructs(species, reactions, c, d)
  # state vector
  x0, id, id2ind = make_species_vector(species, reactions)

  # reactions
  if d <= 8
    r = DenseReactionSystem(reactions, id2ind, c, d)
  else
    r = SparseReactionSystem(reactions, id2ind, c, d)
  end

  # output
  output = initialize_history(c, epochs, trials, time, id2ind)

  return x0, rxn, output
end

function simulate(model::Network, algorithm::Algorithm; sampling_interval::AbstractFloat=1.0, nrlz::Integer=1)

  c = n_species(model)
  d = n_reactions(model)

  species   = species_list(model)
  reactions = reaction_list(model)

  X0, id, id2ind = make_species_vector(species)

  if d <= 8
    r = DenseReactionSystem(reactions, id2ind, c, d)
  else
    r = SparseReactionSystem(reactions, id2ind, c, d)
  end

  t  = end_time(algorithm)

  Xt = deepcopy(X0)

  npts = round(Int, t / sampling_interval + 1)

  if nworkers() == 1
    output = PartialHistory(DenseArray, length(Xt), npts, nrlz, 0.0, t, id2ind)
    serial_simulate(output, Xt, algorithm, X0, r, nrlz)
  else
    output = PartialHistory(SharedArray, length(Xt), npts, nrlz, 0.0, t, id2ind)
    parallel_simulate(output, Xt, algorithm, X0, r)
  end

  return output
end

function serial_simulate(
  output    :: PartialHistory,
  Xt        :: Vector{Int},
  algorithm :: Algorithm,
  X0        :: Vector{Int},
  r         :: AbstractReactionSystem,
  nrlz      :: Integer) # nrlz is encoded in PartialHistory; refactor

  a = propensities(r)

  init!(algorithm, Xt, r)
  for i in 1:nrlz
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
