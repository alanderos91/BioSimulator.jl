mutable struct HybridTauLeapSimulator{M1,M2} <: AbstractSimulator
  exact::M1
  tauleap::M2
  is_critical::Bool
  t::Float64
end

##### constructors #####
function HybridTauLeapSimulator(exact::M1, tauleap::M2) where {M1,M2}
  return HybridTauLeapSimulator(exact, tauleap, false, 0.0)
end

##### main simulation routines #####

# TODO: this is an ugly hack
@inline cumulative_intensity(simulator::HybridTauLeapSimulator) = first(simulator.exact.algorithm.total_rate)

jump_rates(simulator::HybridTauLeapSimulator) = jump_rates(simulator.exact)

@inline function initialize!(simulator::HybridTauLeapSimulator, state, model, tfinal)
  # reset the current system time
  simulator.t = 0.0

  # initialize each simulator
  initialize!(simulator.exact, state, model, tfinal)
  initialize!(simulator.tauleap, state, model, tfinal)

  # check for critical reaction channels
  simulator.is_critical = check_critical(state, model)

  # seed the simulator with events
  # TODO: right now, this seeds both a leap update and jump update
  generate_next_step!(simulator, state)

  return nothing
end

# carry out the next leap
@inline function step!(simulator::HybridTauLeapSimulator, state, model)
  # unpack information
  exact = simulator.exact
  tauleap = simulator.tauleap

  # TODO: right now, each simulator will update by default
  # we should find a way to trigger the correct update
  # based on whether any reaction channel is critical
  if simulator.is_critical
    # perform an update using an exact method
    no_update_step!(exact, state, model)

    simulator.t = exact.t
    tauleap.t = exact.t
  else
    # otherwise it is safe to leap
    no_update_step!(tauleap, state, model)

    simulator.t = tauleap.t
    exact.t = tauleap.t
  end

  # check if simulation is critical
  simulator.is_critical = check_critical(state, model)

  # update the correct simulator
  update!(simulator, state, model)

  # figure out what type of step should be taken next
  generate_next_step!(simulator, state)

  return nothing
end

@inline function update!(simulator::HybridTauLeapSimulator, state, model)
  exact = simulator.exact
  tauleap = simulator.tauleap

  if simulator.is_critical
    update!(exact.algorithm, state, model, 0)
    tauleap.algorithm.total_rate = cumulative_intensity(exact)
  else
    update!(tauleap.algorithm, state, model)
    exact.algorithm.total_rate = cumulative_intensity(tauleap)
  end

  return nothing
end

@inline function generate_next_step!(simulator::HybridTauLeapSimulator, state)
  if simulator.is_critical
    generate_next_jump!(simulator.exact)
  else
    generate_next_leap!(simulator.tauleap, state)
  end

  return nothing
end

@inline function get_new_time(simulator::HybridTauLeapSimulator)
  if simulator.is_critical
    t = get_new_time(simulator.exact)
  else
    t = get_new_time(simulator.tauleap)
  end
  return t
end

# TODO
@inline function check_critical(state, model)
  return any(x -> x < 10, state)
end
