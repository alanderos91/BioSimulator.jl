mutable struct TauLeapSimulator{M,F} <: AbstractSimulator
  algorithm::M
  next_leap_jumps::Vector{Int}
  next_leap_time::Float64
  t::Float64
  execute_leap!::F
end

###### constructors #####
function TauLeapSimulator(algorithm::M, number_jumps) where M
  next_leap_jumps = zeros(Int, number_jumps)

  return TauLeapSimulator{M}(algorithm, next_leap_jumps, 0.0, 0.0)
end

##### main simulation routines #####

# start a new simulation and initialize the stepper with the next event
@inline function initialize!(simulator::TauLeapSimulator, state, model)
  # reset the current system time
  simulator.t = 0.0

  # set the number of jumps to zero
  fill!(simulator.next_leap_jumps, 0)

  # initialize the simulation algorithm
  initialize!(simulator.algorithm, state, model)

  # seed the stepper with the first leap
  generate_next_leap!(simulator)

  return nothing
end

# carry out the next leap
@inline function step!(simulator::TauLeapSimulator, state, model)
  # unpack information
  n = simulator.next_leap_jumps
  s = simulator.next_leap_time

  # update the system time
  simulator.t = get_new_time(simulator)

  # update the system state according to the proposed leap
  simulator.execute_leap!(state, n)

  # update statistics
  # TODO

  # update the simulation algorithm given the leap event
  update!(simulator.algorithm, state, model)

  # generate the next leap
  generate_next_leap!(simulator)

  return nothing
end

@inline function generate_next_leap!(simulator::TauLeapSimulator)
  # unpack
  v = simulator.next_leap_jumps

  # ask the algorithm for the next random leap
  v, s = generate_leap!(v, simulator.algorithm)

  # update the stepper's fields
  simulator.next_leap_jumps = v
  simulator.next_leap_time = s

  return nothing
end

@inline get_new_time(simulator::TauLeapSimulator) = get_new_time(simulator.algorithm, simulator.t, simulator.next_leap_time)
