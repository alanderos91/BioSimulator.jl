abstract Algorithm

steps(x::Algorithm) = x.steps
avg_nsteps(x::Algorithm) = x.avg_nsteps
avg_step_size(x::Algorithm) = x.avg_step_size
done(x::Algorithm, m::Model) = x.t >= x.T || intensity(reactions(m)) == 0.0
time(x::Algorithm) = x.t
end_time(x::Algorithm) = x.T

const DEFAULT_TAGS = [:avg_nsteps, :avg_step_size]

function cumavg(avg, x, n)
    avg = (x + n * avg) / (n + 1)
end

# for use within Algorithm(...)
function compute_statistics!(x::Algorithm, τ::Float64)
    setfield!(x, :avg_step_size, cumavg(avg_step_size(x), τ, steps(x)))
end

# for use at the end of a realization
function compute_statistics!(x::Algorithm, i::Integer)
    setfield!(x, :avg_nsteps, cumavg(avg_nsteps(x), steps(x), i))
end

get_tags(algorithm) = algorithm.tags
