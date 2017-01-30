##### Algorithm interface #####
"""
An abstract type that species a common interface for all algorithms implemented in BioSimulator.

All concrete `Algorithm` types should be immutable and avoid subtyping `Algorithm` directly; see `ExactMethod` and `TauLeapMethod`. The following are normally implemented for `Algorithm`:

- `get_time` : Return the current simulation time.
- `end_time` : Return the termination time specified by the user.
- `done` : Has the simulation terminated? Defaults to `time >= end_time`.
- `init!` : Initialize an `Algorithm` prior to simulation. Memory allocation should happen here because it occurs outside loops in the main simulation routine.
- `step!` : Carry out a simulation step. This should (1) update simulation time, (2) update the species counts, (3) update propensities, and (4) update any relevant data structures.
"""
abstract Algorithm

##### accessors #####
get_time(x::Algorithm) = x.t
end_time(x::Algorithm) = x.end_time

##### termination criterion #####
done(x::Algorithm) = (x.t >= x.end_time)

##### setup outside iteration loop #####
init!(x::Algorithm, Xt, r) = nothing;

##### step through the algorithm #####
step!(x::Algorithm, Xt, r) = nothing;



##### ExactMethod interface #####
"""
An `ExactMethod` is an abstract type for algorithms that are equivalent to Gillespie's method for stochastic simulation.

It subtypes the `Algorithm` type and provides:

- `reset!` : This method defaults to resetting the interal simulation time for an `ExactMethod` after each trial.
"""
abstract ExactMethod <: Algorithm

##### setup inside iteration loop #####
function reset!(algorithm::ExactMethod, a::PVec)
    algorithm.t = 0.0

    return nothing
end

##### TauLeapMethod interface #####
"""
A `TauLeapMethod` is an abstract type for algorithms that utilize the `τ-leaping` strategy to accelerate simulation.

It subtypes the `Algorithm` type and provides:

- `events` : Returns the number of reaction events.
- `reset!` : This method defaults to resetting the interal simulation time for a `TauLeapMethod` after each trial.
"""
abstract TauLeapMethod <: Algorithm

##### accessors #####
events(x::TauLeapMethod) = x.events

##### setup inside iteration loop #####
function reset!(algorithm::TauLeapMethod, a::PVec)
    algorithm.t = 0.0

    return nothing
end
