"""
```
simulate(m::Newtork; with=:ssa, tf=1.0, output=Uniform(), dt=1.0, itr=1, kwargs...)
```

### Arguments
- `m`: The network object

### Optional Arguments
- `with`: The algorithm used in the simulation.
- `tf`: The stopping time for the simulation.
- `output`: The type of output returned by this routine. See [Explicit, Uniform, Histogram].
- `dt`: The step size between updates, if used.
- `itr`: The number of realizations.
- `kwargs`: Optional keyword arguments specific to a particular algorithm. Consult docs for an algorithm for more details.
"""
function simulate(network::Network; with::Symbol=:ssa, T=1.0, output=Uniform(), dt=1.0, itr=1, track=Symbol[], kwargs...)

  if with == :ssa
    algorithm = ssa(T; kwargs...)
  elseif with == :odm
    algorithm = odm(T; kwargs...)
  elseif with == :frm
    algorithm = frm(T; kwargs...)
  elseif with == :nrm
    algorithm = nrm(T; kwargs...)
  elseif with == :sal
    algorithm = sal(T; kwargs...)
  else
    error("$with is an unrecognized algorithm.")
  end

  # Create the main data structure used in the simulation
  Xt, species_id, id2ind = make_species_arr(network.species)
  rv, reaction_id        = make_reaction_arr(network.reactions, id2ind)

  m = Model(network.id, Xt, rv, network.parameters, deepcopy(Xt))

  # Identify which objects to track
  if isempty(track)
      tracked_species = [ i for i in 1:length(species_id) ]
      tracked_reactions = Int[]
  else
      tracked_species   = findin(species_id, track)
      tracked_reactions = findin(reaction_id, track)
  end

  # Create the respective observers
  n = round(Int, T / dt) + 1
  tracker = make_observers(output,
                            species_id,
                            tracked_species,
                            reaction_id,
                            tracked_reactions,
                            Xt,
                            rv,
                            n,
                            itr)

  manager = init_updater(output, tracker, dt, n, itr) # Initialize update manager

  simulate(m, algorithm, output, itr, tracker, manager)
end

function simulate(m::Model, algorithm::Algorithm, output::OutputType, itr, tracker, manager)
  initialize!(algorithm, m)

  for i = 1:itr
    reset!(algorithm)
    reset!(m)
    while !done(algorithm, m)
      update!(output, manager, time(algorithm))   # Record current state
      algorithm(m)
    end
    compute_statistics!(algorithm, i)
    final_update!(output, manager, time(algorithm)) # Record final state
  end

  sdata, pdata = compile_data(tracker)
  mdata = compile_metadata(algorithm, blocksize(manager), itr)
  return SimulationOutput(sdata, pdata, mdata)
end
