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
function simulate(model::Network, algname::T, output_type::Type{Val{O}}=Val{:fixed};
  time::Float64=1.0,
  epochs::Int=1,
  trials::Int=1,
  track_stats::Bool=false,
  kwargs...) where {T,O}

  # build algorithm
  algorithm = build_algorithm(algname, time, track_stats; kwargs...)

  # extract model information
  c = n_species(model)
  d = n_reactions(model)

  species   = species_list(model)
  reactions = reaction_list(model)

  # create simulation data structures
  x0, rxn, id, id2ind = make_datastructs(species, reactions, c, d)

  # get output type
  output = build_output(output_type, c, epochs, trials, time)

  # initialize
  xt = copy(x0)
  init!(algorithm, xt, rxn)

  # delegate trials
  simulate_wrapper!(output, xt, x0, algorithm, rxn)

  result = SimulationSummary(model, algname, algorithm, time, epochs, trials, id2ind, output; kwargs...)

  return result
end

function build_output(::Type{Val{:fixed}}, nspecies, epochs, ntrials, tfinal)
  n = epochs + 1
  tdata = collect(linspace(0.0, tfinal, n))
  output = RegularEnsemble(ntrials, nspecies, epochs)
  for i in eachindex(output)
    @inbounds copy!(output[i].tdata, tdata)
  end
  return output
end

function build_output(::Type{Val{:full}}, nspecies, epochs, ntrials, tfinal)
  return Ensemble(ntrials)
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

function simulate_wrapper!(output, xt, x0, alg, rxn)
  N = Threads.nthreads()
  for i in 1:N
    len = div(length(output), N)
    domain = ((i-1)*len+1):i*len
    simulate_chunk!(output, xt, x0, alg, rxn, domain)
  end
end

function simulate_chunk!(output, Xt, X0, algorithm, reactions, trial_set)
  a = propensities(reactions)
  for trial in trial_set
    copy!(Xt, X0)
    update_all_propensities!(a, reactions, Xt)
    reset!(algorithm, a)

    xw = output[trial]

    simulate!(xw, Xt, algorithm, reactions)
  end

  return output
end

function simulate!(xw :: RegularPath, Xt, algorithm, reactions)
  epoch = 1
  while !done(algorithm)
    epoch = update!(xw, algorithm.t, Xt, epoch)
    step!(algorithm, Xt, reactions)
  end
  update!(xw, algorithm.t, Xt, epoch)

  return xw
end

function simulate!(xw :: SamplePath, Xt, algorithm, reactions)
  while !done(algorithm)
    update!(xw, algorithm.t, Xt)
    step!(algorithm, Xt, reactions)
  end

  return xw
end
