"""
```
simulate(m::Newtork; method=:SSA, time=1.0, output=:explicit, sampling_interval=1.0, realizations=1, kwargs...)
```

### Arguments
- `m`: The network object

### Optional Arguments
- `method`: The algorithm used in the simulation.
- `time`: The stopping time for the simulation.
- `output`: The type of output to record. The `:explicit` option will record simulation history after every algorithm step, whereas `:fixed` will record the history after according to `sampling_interval`.
- `sampling_interval`: The time between updates when `output=:fixed`. The default is to record after every unit of time, according to the scale implicit in a given model `m`.
- `realizations`: The number of Monte Carlo realizations.
- `kwargs`: Optional keyword arguments specific to a particular algorithm. Refer to the appropriate algorithm for details.
"""
function simulate(network::Network; method::Symbol=:SSA, time=1.0, output=:fixed, sampling_interval=1.0, realizations=1, track=Symbol[], kwargs...)

  if method == :SSA
    algorithm = ssa(time; kwargs...)
  elseif method == :ODM
    algorithm = odm(time; kwargs...)
  elseif method == :FRM
    algorithm = frm(time; kwargs...)
  elseif method == :NRM
    algorithm = nrm(time; kwargs...)
  elseif method == :SAL
    algorithm = sal(time; kwargs...)
  else
    error("$method is an unrecognized algorithm.")
  end

  otype = (output == :explicit) ? Explicit() : Uniform()

  # Create the main data structure used in the simulation
  Xt, species_id, id2ind = make_species_arr(species_list(network))
  rs = reaction_system(reaction_list(network), id2ind)
  reaction_id = Symbol[]

  m = Model(network.id, Xt, rs, parameter_list(network), deepcopy(Xt))

  # Identify which objects to track
  if isempty(track)
      tracked_species = [ i for i in 1:length(species_id) ]
      tracked_reactions = Int[]
  else
      tracked_species   = findin(species_id, track)
      tracked_reactions = findin(reaction_id, track)
  end

  Xt_history = initialize_history(otype, Xt, id2ind, realizations, Int(time / sampling_interval) + 1)

  simulate(m, algorithm, otype, realizations, Xt_history)
end

function simulate(m::Model, algorithm::Algorithm, output::OutputType, realizations, Xt_history)
  initialize!(algorithm, m)

  Xt = m.Xt
  rs = m.rs
  p  = m.parameters
  X₀ = m.X₀

  for i = 1:realizations
    reset!(m)
    reset!(algorithm, m)
    while time(algorithm) < end_time(algorithm)
      #@code_warntype update!(Xt_history, i, time(algorithm), Xt)
      #update!(Xt_history, i, time(algorithm), Xt)   # Record current state
      step!(algorithm, Xt, rs, p)
    end
    compute_statistics!(algorithm, i)
    #update!(Xt_history, i, time(algorithm), Xt) # Record final state
  end

  return Xt_history
end

function initialize_history(::Explicit, Xt, id2ind, nrlz, nstates)
  return SimData(Xt, nrlz, id2ind)
end

function initialize_history(::Uniform, Xt, id2ind, nrlz, nstates)
  return SimData(Xt, nrlz, id2ind, nstates)
end
