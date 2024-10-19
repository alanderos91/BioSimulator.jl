mutable struct TauLeapSimulator{M,F} <: AbstractSimulator
  algorithm::M
  next_leap_jumps::Vector{Int}
  next_leap_time::Float64
  t::Float64
  tfinal::Float64
  execute_leap!::F
end

###### constructors #####
function TauLeapSimulator(algorithm::M, number_jumps, execute_leap!::F) where {M,F}
  next_leap_jumps = zeros(Int, number_jumps)

  return TauLeapSimulator{M,F}(algorithm, next_leap_jumps, 0.0, 0.0, 0.0, execute_leap!)
end

##### main simulation routines #####

@inline cumulative_intensity(simulator::TauLeapSimulator) = first(simulator.algorithm.total_rate)

jump_rates(simulator::TauLeapSimulator) = simulator.algorithm.rates

# start a new simulation and initialize the stepper with the next event
@inline function initialize!(simulator::TauLeapSimulator, state, model, tfinal)
  # reset the current system time
  simulator.t = 0.0

  simulator.tfinal = tfinal

  # set the number of jumps to zero
  fill!(simulator.next_leap_jumps, 0)

  # initialize the simulation algorithm
  initialize!(simulator.algorithm, state, model)

  # seed the stepper with the first leap
  generate_next_leap!(simulator, state)

  return nothing
end

# carry out the next leap
@inline function step!(simulator::TauLeapSimulator, state, model)
  # unpack information
  jumps = simulator.next_leap_jumps

  # update the system time
  simulator.t = get_new_time(simulator)

  # update the system state according to the proposed leap
  # NOTE: at this point it is assumed that the leap is valid
  simulator.execute_leap!(state, jumps)

  # update statistics
  # TODO

  # update the simulation algorithm given the leap event
  update!(simulator.algorithm, state, model)

  # generate the next leap
  generate_next_leap!(simulator, state)

  return nothing
end

@inline function generate_next_leap!(simulator::TauLeapSimulator, state)
  # set a limit on how large a leap may be
  Δt_max = simulator.tfinal - simulator.t

  # ask the algorithm for the next random leap
  v, s = generate_leap!(simulator.algorithm, simulator.next_leap_jumps, state, Δt_max)

  # update the stepper's fields
  copyto!(simulator.next_leap_jumps, v)
  simulator.next_leap_time = s

  return nothing
end

@inline get_new_time(simulator::TauLeapSimulator) = get_new_time(simulator.algorithm, simulator.t, simulator.next_leap_time)

##### helper functions #####

@inline function no_update_step!(simulator::TauLeapSimulator, state, model)
  # update the system time
  simulator.t = get_new_time(simulator)

  # update the system state according to the proposed leap
  # NOTE: at this point it is assumed that the leap is valid
  simulator.execute_leap!(state, simulator.next_leap_jumps)

  return nothing
end
