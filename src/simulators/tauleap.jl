mutable struct TauLeapSimulator{M} <: AbstractSimulator
  algorithm::M
  next_leap_jumps::Vector{Int}
  next_leap_time::Float64
  t::Float64
end

###### constructors #####
function TauLeapSimulator(algorithm::M) where M
  return TauLeapSimulator{M}(algorithm, Int[], 0.0, 0.0)
end

##### main simulation routines #####

# start a new simulation and initialize the stepper with the next event
@inline function initialize!(simulator::TauLeapSimulator, state, model)
  # reset the current system time
  simulator.t = 0.0

  # initialize the simulation algorithm
  initialize!(simulator.algorithm, state, model)

  # seed the stepper with the first leap
  generate_next_leap!(simulator)

  return nothing
end

# carry out the next leap
@inline function step!(simulator::TauLeapSimulator, staet, model)
  # unpack information
  v = simulator.next_leap_jumps
  s = simulator.next_leap_time

  # update the system time
  simulator.t = get_new_time(simulator)

  # update the system state according to the proposed leap
  @inbounds @fastmath execute_leap!(state, model, v)

  # update statistics
  # TODO

  # update the simulation algorithm given the leap event
  update!(simulator.algorithm, state, model, v)

  # generate the next leap
  generate_next_leap!(simulator)

  return nothing
end

@inline function generate_next_leap!(simulator::TauLeapSimulator)
  # unpack
  v = simulator.next_leap_jumps
  
  # ask the algorithm for the next random leap
  v, s = generate_leap!(simulator.algorithm)

  # update the stepper's fields
  simulator.next_leap_jumps = v
  simulator.next_leap_time = s

  return nothing
end

@inline get_new_time(simulator::TauLeapSimulator) = get_new_time(simulator.algorithm, simulator.t, simulator.next_leap_time)