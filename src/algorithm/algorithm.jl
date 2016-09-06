##### Algorithm interface #####

abstract Algorithm

const DEFAULT_TIME = 1.0

##### accessors #####
get_time(x::Algorithm) = x.t
end_time(x::Algorithm) = x.end_time
get_tags(x::Algorithm) = x.tags

##### termination criterion #####
done(x::Algorithm) = x.t >= x.end_time

##### setup outside iteration loop #####
init!(x::Algorithm, Xt, r) = nothing;

##### step through the algorithm #####
step!(x::Algorithm, Xt, r) = nothing;



##### ExactMethod interface #####

abstract ExactMethod <: Algorithm

const DEFAULT_EXACT = [:avg_nsteps, :avg_stepsz]

##### accessors #####
steps(x::ExactMethod)      = x.nsteps
avg_nsteps(x::ExactMethod) = x.avg_nsteps
avg_stepsz(x::ExactMethod) = x.avg_stepsz

function nsteps!(algorithm::ExactMethod)
    algorithm.nsteps += 1
    return nothing
end

##### setup inside iteration loop #####
function reset!(x::ExactMethod, a::PVec)
    setfield!(x, :t, 0.0)
    setfield!(x, :nsteps, 0)
    return;
end

function compute_statistics!(x::ExactMethod, τ::AbstractFloat)
    setfield!(x, :avg_stepsz, cumavg(avg_stepsz(x), τ, steps(x)))

    return nothing
end

function compute_statistics!(x::ExactMethod, i::Integer)
    setfield!(x, :avg_nsteps, cumavg(avg_nsteps(x), steps(x), i - 1))

    return nothing
end

##### TauLeapMethod interface #####

abstract TauLeapMethod <: Algorithm

const DEFAULT_TAULEAP = [:avg_nsteps,
                         :avg_stepsz,
                         :avg_nssa,
                         :avg_nleaps,
                         :avg_ssa_step,
                         :avg_leap_step,
                         :avg_neg_excursions]

##### accessors #####
steps(x::TauLeapMethod)              = x.ssa_steps + x.leap_steps
ssa_steps(x::TauLeapMethod)          = x.ssa_steps
leap_steps(x::TauLeapMethod)         = x.leap_steps
neg_excursions(x::TauLeapMethod)     = x.neg_excursions

avg_nsteps(x::TauLeapMethod)         = x.avg_nsteps
avg_stepsz(x::TauLeapMethod)         = x.avg_stepsz
avg_nssa(x::TauLeapMethod)           = x.avg_nssa
avg_nleaps(x::TauLeapMethod)         = x.avg_nleaps
avg_ssa_step(x::TauLeapMethod)       = x.avg_ssa_step
avg_leap_step(x::TauLeapMethod)      = x.avg_leap_step
avg_neg_excursions(x::TauLeapMethod) = x.avg_neg_excursions

events(x::TauLeapMethod) = x.events

function nsteps!(algorithm::TauLeapMethod, isleap::Bool)
    if isleap
        algorithm.leap_steps += 1
    else
        algorithm.ssa_steps += 1
    end

    return nothing
end

##### setup inside iteration loop #####
function reset!(x::TauLeapMethod, a::PVec)
    setfield!(x, :t, 0.0)
    setfield!(x, :ssa_steps, 0)
    setfield!(x, :leap_steps, 0)
    setfield!(x, :neg_excursions, 0)

    return nothing
end

function compute_statistics!(x::TauLeapMethod, τ::AbstractFloat, isleap::Bool)
    setfield!(x, :avg_stepsz, cumavg(avg_stepsz(x), τ, steps(x)))

    if isleap
        setfield!(x, :avg_leap_step, cumavg(avg_leap_step(x), τ, leap_steps(x)))
    else
        setfield!(x, :avg_ssa_step, cumavg(avg_ssa_step(x), τ, ssa_steps(x)))
    end

    return nothing
end

function compute_statistics!(x::TauLeapMethod, i::Integer)
    setfield!(x, :avg_nsteps, cumavg(avg_nsteps(x), steps(x), i - 1))
    setfield!(x, :avg_nssa, cumavg(avg_nssa(x), ssa_steps(x), i - 1))
    setfield!(x, :avg_nleaps, cumavg(avg_nleaps(x), leap_steps(x), i - 1))

    return nothing
end

##### misc #####

# cumulative average
cumavg(avg, x, n::Integer) = (x + n * avg) / (n + 1)
