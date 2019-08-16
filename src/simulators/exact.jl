mutable struct ExactSimulator{M} <: AbstractSimulator
  algorithm::M
  next_jump_index::Int
  next_jump_time::Float64
  t::Float64
end

##### constructors ######
function ExactSimulator(algorithm::M) where M
  return ExactSimulator{M}(algorithm, 0, 0.0, 0.0)
end

##### main simulation routines #####

@inline cumulative_intensity(simulator::ExactSimulator) = first(simulator.algorithm.total_rate)

# start a new simulation and initialize the stepper with the next event
@inline function initialize!(simulator::ExactSimulator, state, model, tfinal)
  # reset the current system time
  simulator.t = 0.0

  # initialize the simulation algorithm
  initialize!(simulator.algorithm, state, model)

  # seed the stepper with the first event
  generate_next_jump!(simulator)

  return nothing
end

# carry out the next jump event
@inline function step!(simulator::ExactSimulator, state, model)
  # unpack information
  j = simulator.next_jump_index
  s = simulator.next_jump_time

  # update the system time
  simulator.t = get_new_time(simulator)

  # update the system state according to jump j
  @inbounds @fastmath execute_jump!(state, model, j)

  # update statistics
  # TODO

  # update the simulation algorithm given that jump j occured
  update!(simulator.algorithm, state, model, j)

  # generate the next jump
  generate_next_jump!(simulator)

  return nothing
end

# generate the next (j,s) pair by invoking the stepper's event sampler
@inline function generate_next_jump!(simulator::ExactSimulator)
  # ask the algorithm for the next random jump
  j, s = generate_jump(simulator.algorithm)

  # update the stepper's fields
  simulator.next_jump_index = j
  simulator.next_jump_time = s

  return nothing
end

@inline get_new_time(simulator::ExactSimulator) = get_new_time(simulator.algorithm, simulator.t, simulator.next_jump_time)

##### helper functions #####

function no_update_step!(simulator::ExactSimulator, state, model)
  # unpack information
  j = simulator.next_jump_index
  s = simulator.next_jump_time

  # update the system time
  simulator.t = get_new_time(simulator)

  # update the system state according to jump j
  @inbounds @fastmath execute_jump!(state, model, j)

  return nothing
end
