"""
```
simulate{T}(model, [algorithm::Type{T}=SSA];
  time   :: AbstractFloat=1.0,
  epochs :: Integer=1,
  trials :: Integer=1,
  algvars...)
```

Simulate a `Network` using the given `algorithm`.

The simulation routine will run until the termination `time` and record the system state at evenly spaced `epochs`. This repeats for the given number of `trials`.

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
  time::Float64=1.0,
  epochs::Int=1,
  trials::Int=1,
  kwargs...)

  # extract model information
  c = n_species(model)
  d = n_reactions(model)

  species = species_list(model)
  reactions = reaction_list(model)

  # create simulation data structures
  x0, rxn, id, id2ind = make_datastructs(species, reactions, c, d)
  t_index = linspace(0.0, time, epochs + 1)

  # initialize
  xt     = copy(x0)
  alg    = algorithm(time; kwargs...)
  output = SimData(SharedArray(eltype(x0), c, epochs + 1, trials))

  init!(alg, xt, rxn)

  # simulation
  @sync for pid in procs(output.data)
    @async remotecall_fetch(pid, simulate_shared_chunk!, output, xt, x0, alg, rxn, t_index)
  end

  return output
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

# generate a single trajectory
function simulate!(
    output    :: SimData,
    Xt        :: Vector{Int},
    algorithm :: Algorithm,
    reactions :: AbstractReactionSystem,
    t_index   :: Range,
    trial     :: Integer
  )
  epoch = 1
  while !done(algorithm)
    t     = get_time(algorithm)
    epoch = update!(output, Xt, t, t_index, epoch, trial)

    step!(algorithm, Xt, reactions)
  end
  t = get_time(algorithm)
  update!(output, Xt, t, t_index, epoch, trial)

  return output
end

# this retrieves trials assigned to process
function partition(output :: SimData)
  q   = get_data(output)
  idx = indexpids(q)
  if idx == 0
    # This worker is not assigned a piece
    return 1:0, 1:0
  end
  nchunks = length(procs(q))
  splits  = [round(Int, s) for s in linspace(0,size(q,3),nchunks+1)]
  splits[idx]+1:splits[idx+1]
end

# iterate over trials assigned to process
function simulate_chunk!(
  output    :: SimData,
  Xt        :: Vector{Int},
  X0        :: Vector{Int},
  algorithm :: Algorithm,
  reactions :: AbstractReactionSystem,
  t_index   :: Range,
  trial_set :: UnitRange
)
  a = propensities(reactions)

  for trial in trial_set
    copy!(Xt, X0)
    update_all_propensities!(a, reactions, Xt)
    reset!(algorithm, a)

    simulate!(output, Xt, algorithm, reactions, t_index, trial)
  end

  return output
end

# the wrapper
@inline function simulate_shared_chunk!(
  output    :: SimData,
  Xt        :: Vector{Int},
  X0        :: Vector{Int},
  algorithm :: Algorithm,
  reactions :: AbstractReactionSystem,
  t_index   :: Range
)
  simulate_chunk!(output, Xt, X0, algorithm, reactions, t_index, partition(output))
end
