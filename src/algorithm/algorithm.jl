##### Algorithm interface #####

abstract Algorithm

##### accessors #####
time(x::Algorithm)     = x.t
end_time(x::Algorithm) = x.T
get_tags(x::Algorithm) = x.tags

##### termination criterion #####
done(x::Algorithm, m::Model) = x.t >= x.T || intensity(reactions(m)) == 0.0

##### setup outside iteration loop #####
initialize!(x::Algorithm, m::Model) = return;

#### setup inside iteration loop #####
reset!(x::Algorithm, m::Model) = return;

##### step through the algorithm #####
step!(x::Algorithm, Xt, rs, p) = return;




##### ExactMethod interface #####

abstract ExactMethod <: Algorithm

const DEFAULT_EXACT = [:avg_nsteps, :avg_step_size]

##### accessors #####
steps(x::ExactMethod)         = x.steps
avg_nsteps(x::ExactMethod)    = x.avg_nsteps
avg_step_size(x::ExactMethod) = x.avg_step_size

##### setup inside iteration loop #####
function reset!(x::ExactMethod, m::Model)
    setfield!(x, :t, 0.0)
    setfield!(x, :steps, 0)
    return;
end

# use within the step! method
function compute_statistics!(x::ExactMethod, τ::Float64)
    setfield!(x, :avg_step_size, cumavg(avg_step_size(x), τ, steps(x)))
    return;
end

# use at the end of a realization
function compute_statistics!(x::ExactMethod, i::Integer)
    setfield!(x, :avg_nsteps, cumavg(avg_nsteps(x), steps(x), i - 1))
    return;
end




##### TauLeapMethod interface #####

abstract TauLeapMethod <: Algorithm

const DEFAULT_TAULEAP = [:avg_nsteps, :avg_step_size, :avg_nssa, :avg_nleaps, :avg_ssa_step, :avg_leap_step, :avg_neg_excursions]

##### accessors #####
steps(x::TauLeapMethod)              = x.ssa_steps + x.leap_steps
ssa_steps(x::TauLeapMethod)          = x.ssa_steps
leap_steps(x::TauLeapMethod)         = x.leap_steps
neg_excursions(x::TauLeapMethod)     = x.neg_excursions

avg_nsteps(x::TauLeapMethod)         = x.avg_nsteps
avg_step_size(x::TauLeapMethod)      = x.avg_step_size
avg_nssa(x::TauLeapMethod)           = x.avg_nssa
avg_nleaps(x::TauLeapMethod)         = x.avg_nleaps
avg_ssa_step(x::TauLeapMethod)       = x.avg_ssa_step
avg_leap_step(x::TauLeapMethod)      = x.avg_leap_step
avg_neg_excursions(x::TauLeapMethod) = x.avg_neg_excursions

events(x::TauLeapMethod) = x.events

##### setup inside iteration loop #####
function reset!(x::TauLeapMethod, m::Model)
    setfield!(x, :t, 0.0)
    setfield!(x, :ssa_steps, 0)
    setfield!(x, :leap_steps, 0)
    setfield!(x, :neg_excursions, 0)
    return;
end

# use within the step! method
function compute_statistics!(x::TauLeapMethod, τ::Float64)
    setfield!(x, :avg_step_size, cumavg(avg_step_size(x), τ, steps(x)))
end

# use at the end of a realization
function compute_statistics!(x::TauLeapMethod, i::Integer)
    setfield!(x, :avg_nsteps, cumavg(avg_nsteps(x), steps(x), i - 1))
    setfiedl!(x, :avg_neg_excursions, cumavg(avg_neg_excursions(x), steps(x), i -1))
end

##### misc #####

# cumulative average
function cumavg(avg, x, n)
    avg = (x + n * avg) / (n + 1)
end
