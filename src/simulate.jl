function parse_model(network::Network)
  state = [ x.population for (key, x) in species_list(network) ]
  model = ReactionSystem(network)

  return state, model
end

function simulate(network::Network, algname::SimulationAlgorithm; tfinal=0.0, rates_cache = HasRates)
  # build the internal representation of our stochastic process
  initial_state, model = parse_model(network)

  # feedforward down the chain...
  return simulate(initial_state, model, algname, tfinal, rates_cache)
end

function simulate(initial_state, model, algname, tfinal, rates_cache)
  # copy state
  state = copy(initial_state)

  # build the simulator
  simulator = build_simulator(algname, state, model, rates_cache)

  # build the output data
  output = build_output(state, model)

  initialize_datastructs!(state, model)

  # feedforward down the chain...
  simulate!(simulator, state, model, tfinal, output)
end

function simulate!(simulator, state, model, tfinal, output)
  initialize!(simulator, state, model)

  while simulator.t < tfinal && first(simulator.algorithm.total_rate) > 0
    tnew = get_new_time(simulator)

    if tnew <= tfinal
      step!(simulator, state, model)
      update!(output, simulator.t, state)
    else
      simulator.t = tfinal
    end
  end
  update!(output, simulator.t, state)

  return output
end

##### temporary hacks

function build_output(state, model)
  xw = SamplePath([copy(state)], [0.0])
  sizehint!(xw, 1_000)

  return xw
end

function build_output(state::Lattice, model)
  xw = ([Configuration(state)], [0.0])
  sizehint!(xw[1], 1_000)
  sizehint!(xw[2], 1_000)

  return xw
end

function update!(xw, t, x::Lattice)
  push!(xw[1], Configuration(x))
  push!(xw[2], t)

  return xw
end

initialize_datastructs!(state, model) = nothing

function initialize_datastructs!(lattice::Lattice, model::InteractingParticleSystem)
  # unpack information
  number_init = number_sites(lattice)
  site = lattice.site
  nbtype = topology(lattice)
  enum = model.enumeration
  class = enum.class
  isactive = enum.isactive
  pair_to_classes = enum.pair_to_classes
  dummy_composition = enum.dummy_composition

  for s in eachindex(class)
    empty!(class[s])
  end

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
