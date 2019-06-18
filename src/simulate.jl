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

  # feedforward down the chain...
  simulate!(simulator, state, model, tfinal, output)
end

function simulate!(simulator, state, model, tfinal, output)
  initialize!(simulator, state, model)

  update!(output, simulator.t, state)

  while simulator.t < tfinal && first(simulator.algorithm.total_rate) > 0
    tnew = get_new_time(simulator)

    if tnew <= tfinal
      step!(simulator, state, model)
    else
      simulator.t = tfinal
    end

    update!(output, simulator.t, state)
  end

  return output
end

# this is temporary
function build_output(state, model)
  xw = SamplePath{Vector{Float64}, Vector{Vector{Int}}}()
  sizehint!(xw, 1_000)

  return xw
end
