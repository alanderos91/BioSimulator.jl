"""
```
parse_model(network::Network)
```

Construct a `state <: Vector{Int}` and `model <: ReactionSystem` from a `Network` object.

Use this in `simulate` to avoid constructing internal data structures over multiple samples.
"""
function parse_model(network::Network)
    state = [ x.population for (key, x) in species_list(network) ]
    model = ReactionSystem(network)

    return state, model
end

##### user-facing interface #####

"""
```
simulate(network::Network, algname::SimulationAlgorithm;
            [tfinal = 0.0],
            [rates_cache = HasRates],
            [save_points = nothing],
            [save_function = save_state],
            [ntrials = nothing])
```

Simulate a `Network` of interacting populations.
Note that simulations may terminate early if the cumulative intensity reaches `0` (that is, reaching an absorbing state).

#### Keyword Arguments

- `tfinal`: The final simulation time.
- `rates_cache`: Indicates the type of information stored in a `rates` cache. If `HasRates` is chosen, store the actual rates. If `HasSums` is chosen, store partial sums of the `rates`. This effectively toggles between linear and binary searches, provided the algorithm supports the option.
- `save_points`: An indexable collection indicating time points to sample and record system state. The default `nothing` forces saving after every reaction event.
- `save_function`: A three argument function that maps the three arguments `(simulator, state, model)` to data recorded in a `SamplePath`. Default behavior is to record each species.
- `ntrials`: Specifies the number of times to simulate the given model. The default `nothing` generates a single sample path.
"""
function simulate(network::Network, algname::SimulationAlgorithm;
    tfinal=0.0,
    rates_cache = HasRates,
    save_points = nothing,
    save_function::funcT = save_state,
    ntrials = nothing
    ) where funcT
    # build the internal representation of our stochastic process
    initial_state, model = parse_model(network)

    # feedforward down the chain...
    return simulate(initial_state, model, algname, tfinal, rates_cache, save_points, save_function, ntrials)
end

"""
```
simulate(state, model, algname::SimulationAlgorithm;
            [tfinal = 0.0],
            [rates_cache = HasRates],
            [save_points = nothing],
            [save_function = save_state],
            [ntrials = nothing])
```

Simulate a `model` with the given initial `state`.
Note that simulations may terminate early if the cumulative intensity reaches `0` (that is, reaching an absorbing state).

#### Arguments

For well-mixed systems:

- `state`: A `Vector{Int}` constructed by `parse_model`.
- `model`: A `ReactionSystem` type constructed by `parse_model`.

For lattice-based systems:

- `state`: A `Lattice` type whose topology comptabile is compatible with `model`.
- `model`: An `InteractingParticleSystem` type built from `@enumerate_with_sclass`.

#### Keyword Arguments

- `tfinal`: The final simulation time.
- `rates_cache`: Indicates the type of information stored in a `rates` cache. If `HasRates` is chosen, store the actual rates. If `HasSums` is chosen, store partial sums of the `rates`. This effectively toggles between linear and binary searches, provided the algorithm supports the option.
- `save_points`: An indexable collection indicating time points to sample and record system state. The default `nothing` forces saving after every reaction event.
- `save_function`: A three argument function that maps the three arguments `(simulator, state, model)` to data recorded in a `SamplePath`. Default behavior is to store population counts for well-mixed systems or `Configuration` objects for lattice-based systems.
- `ntrials`: Specifies the number of times to simulate the given model. The default `nothing` generates a single sample path.
"""
function simulate(initial_state, model, algname::SimulationAlgorithm;
    tfinal = 0.0,
    rates_cache = HasRates,
    save_points = nothing,
    save_function::funcT = save_state,
    ntrials = nothing
    ) where funcT
    # feedforward down the chain...
    return simulate(initial_state, model, algname, tfinal, rates_cache, save_points, save_function, ntrials)
end

##### internals #####

function simulate(initial_state, model, algname, tfinal, rates_cache, save_points, save_function::F, ntrials) where F
    # copy state
    state = copy(initial_state)

    # build the simulator
    simulator = build_simulator(algname, state, model, rates_cache)
    # TODO: figure out a way of initializing so that output is built correctly
    initialize!(simulator, state, model, tfinal)

    # build the output data
    output = build_output(save_function, simulator, state, model, ntrials)

    # feedforward down the chain...
    simulate!(simulator, state, initial_state, model, tfinal, output, save_points, save_function, ntrials)
end

# case: single trajectory
function simulate!(simulator, state, initial_state, model, tfinal, output, save_points, save_function::F, ::Nothing) where F
    initialize_datastructs!(state, initial_state, model)
    simulate!(simulator, state, model, tfinal, output, save_points, save_function)

    return output
end

# case: ensemble simulation
function simulate!(simulator, state, initial_state, model, tfinal, output, save_points, save_function::F, ntrials::Integer) where F
    for k in 1:ntrials
        initialize_datastructs!(state, initial_state, model)
        sample_path = output[k]
        simulate!(simulator, state, model, tfinal, sample_path, save_points, save_function)
    end

    return output
end

function simulate!(simulator, state, model, tfinal, output, save_points, save_function::F) where F
    initialize!(simulator, state, model, tfinal)

    in_timespan = simulator.t < tfinal
    transient = cumulative_intensity(simulator) > 0
    running = in_timespan && transient

    while running
        tnew = get_new_time(simulator)

        if tnew < tfinal
            step!(simulator, state, model)
        else
            simulator.t = tfinal
        end
        # update output
        t = simulator.t
        update_samplepath!(save_function, output, simulator, t, state, model, save_points)

        # update flags
        in_timespan = simulator.t < tfinal
        transient = cumulative_intensity(simulator) > 0
        running = in_timespan && transient
    end

    # add end point in case of extinction
    if !transient
        update_samplepath!(save_function, output, simulator, tfinal, state, model, save_points)
    end

    return output
end

##### temporary hacks

initialize_datastructs!(state, initial_state, model) = copyto!(state, initial_state)

function initialize_datastructs!(lattice::Lattice, initial_lattice::Lattice, model::InteractingParticleSystem)
    # copy fields for lattice...
    empty!(lattice.site)
    for s in initial_lattice.site
        new_site = Site(s.label, State(s.state.ptype), tuple(s.coord...))
        push!(lattice.site, new_site)
    end

    empty!(lattice.coord_order)
    foreach(s -> push!(lattice.coord_order, s), lattice.site)
    sort!(lattice.coord_order,
    alg = Base.Sort.DEFAULT_STABLE,
    by = coordinates,
    lt = <)

    foreach(nb -> empty!(nb), lattice.neighbors)

    # unpack information
    number_init = number_sites(lattice)
    site = lattice.site
    nbtype = topology(lattice)
    enum = model.enumeration
    class = enum.class
    isactive = enum.isactive
    pair_to_classes = enum.pair_to_classes
    dummy_composition = enum.dummy_composition

    # make sure sample classes are empty
    foreach(c -> empty!(c), class)

    # add missing open sites
    for i in 1:number_init
        x = site[i]
        for coord in eachdir(nbtype, x)
            if !istracked(lattice, coord)
                y = spawn_new_site(lattice, coord)
                add_site!(lattice, y)
            end
        end
    end

    # build neighborhoods
    build_neighborhoods!(nbtype, lattice)

    # assign sites to classes in the enumeration
    for x in site
        # build the neighborhood composition around x
        build_composition!(dummy_composition, lattice, x)

        # what sample classes does x belong to?
        l = get_ptype(x)
        k = get_nbclass_index(enum, dummy_composition)

        if isactive[l]
            classes = pair_to_classes[(l,k)]

            # add x to each assigned class
            for s in classes
                add_class_member!(class[s], x)
            end
        end

        # change the neighborhood class of x
        change_neighbor_class!(x, k)
    end

    # clean up
    fill!(dummy_composition, 0)
    sort!(lattice.coord_order,
    alg = Base.Sort.DEFAULT_STABLE,
    by = coordinates,
    lt = <)

    return nothing
end
