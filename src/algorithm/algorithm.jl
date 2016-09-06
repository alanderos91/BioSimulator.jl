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

##### setup inside iteration loop #####
function reset!(x::ExactMethod, a::PVec)
    setfield!(x, :t, 0.0)
    setfield!(x, :nsteps, 0)
    return;
end

function update_statistics!(x::ExactMethod, τ::AbstractFloat)
    x.nsteps += 1
    fit!(avg_stepsz(x), τ)

    return nothing
end

function update_statistics!(x::ExactMethod)
    fit!(avg_nsteps(x), steps(x))

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

##### setup inside iteration loop #####
function reset!(x::TauLeapMethod, a::PVec)
    setfield!(x, :t, 0.0)
    setfield!(x, :ssa_steps, 0)
    setfield!(x, :leap_steps, 0)
    setfield!(x, :neg_excursions, 0)

    return nothing
end

function update_statistics!(x::TauLeapMethod, τ::AbstractFloat, isleap::Bool)
    fit!(avg_stepsz(x), τ)

    if isleap
        x.leap_steps += 1
        fit!(avg_leap_step(x), τ)
    else
        x.ssa_steps += 1
        fit!(avg_ssa_step(x), τ)
    end

    return nothing
end

function update_statistics!(x::TauLeapMethod)
    fit!(avg_nsteps(x), steps(x))
    fit!(avg_nssa(x), ssa_steps(x))
    fit!(avg_nleaps(x), leap_steps(x))
    fit!(avg_neg_excursions(x), neg_excursions(x))

    return nothing
end

##### misc #####

# cumulative average
cumavg(avg, x, n::Integer) = (x + n * avg) / (n + 1)
