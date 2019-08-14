mutable struct HybridTauLeapSimulator{M1,M2} <: AbstractSimulator
  exact_algorithm::M1
  next_jump_index::Int
  next_jump_time::Float64
  tauleap_algorithm::M2
  next_leap_jumps::Vector{Int}
  next_leap_time::Float64
  is_critical::Bool
  t::Float64
  tfinal::Float64
  execute_leap!::F
end

##### constructors #####
function HybridTauLeapSimulator(exact::M1, tauleap::M2, execute_leap!) where {M1,M2}
  return HybridTauLeapSimulator(exact, 0, 0.0, tauleap, ?, 0.0, false, 0.0, 0.0, execute_leap!)
end

##### main simulation routines #####

@inline function initialize!(simulator::HybridTauLeapSimulator, state, model, tfinal)
  # reset the current system time
  simulator.t = 0.0

  # initialize each algorithm
  initialize!(simulator.exact_algorithm, state, model, tfinal)
  initialize!(simulator.tauleap_algorithm, state, model, tfinal)

  # seed the simulator with events

  return nothing
end

# carry out the next leap
@inline function step!(simulator::HybridTauLeapSimulator, state, model)
  # unpack information
  exact_algorithm = simulator.exact_algorithm
  tauleap_algorithm = simulator.tauleap_algorithm
  execute_leap! = simulator.execute_leap!

  if simulator.is_critical
    # perform an update using an exact method
    simulator.t = get_new_time(exact_algorithm, simulator.t, simulator.next_jump_time)
    execute_jump!(state, model, simulator.next_jump_index)
  else
    # otherwise it is safe to leap
    simulator.t = get_new_time(tauleap_algorithm, simulator.t, simulator.next_leap_time)
    execute_leap!(state, simulator.next_leap_jumps)
  end

  # update the simulator; need to determine if system is critical
  update!(simulator.algorithm, state, model)

  # figure out what type of step should be taken next
  generate_next_leap!(simulator)

  return nothing
end

@inline function generate_next_leap!(simulator::HybridTauLeapSimulator)

end

@inline get_new_time(simulator::HybridTauLeapSimulator)

end

