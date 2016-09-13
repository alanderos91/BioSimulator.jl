##### Algorithm interface #####

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

abstract ExactMethod <: Algorithm

##### setup inside iteration loop #####
function reset!(algorithm::ExactMethod, a::PVec)
    algorithm.t = 0.0

    return nothing
end

##### TauLeapMethod interface #####

abstract TauLeapMethod <: Algorithm

##### accessors #####
events(x::TauLeapMethod) = x.events

##### setup inside iteration loop #####
function reset!(algorithm::TauLeapMethod, a::PVec)
    algorithm.t = 0.0

    return nothing
end
