type FRM <: Algorithm
    # parameters
    T::Float64

    # state variables
    t::Float64
    steps::Int

    # statistics
    avg_nsteps::Float64
    avg_step_size::Float64

    # metadata tags
    tags::Vector{Symbol}

    function FRM(T)
        new(T, 0.0, 0, 0.0, 0.0, DEFAULT_TAGS)
    end
end

function frm(T; na...)
    return FRM(T)
end

initialize!(x::FRM, m::Model) = return;

function reset!(x::FRM, m::Model)
    setfield!(x, :t, 0.0)
    setfield!(x, :steps, 0)
    return;
end

function call(x::FRM, m::Model)
    compute_propensities!(m)
    τ = frm_update!(m, x)

    # update algorithm variables
    setfield!(x, :t,     time(x) + τ)
    setfield!(x, :steps, steps(x) + 1)

    compute_statistics!(x, τ)
end

function frm_update!(m::Model, x)
    t = x.t
    T = x.T

    rxns = reactions(m)

    τ = Inf; μ = 0
    for j in eachindex(rxns)
        τj = rand(Exponential(1/rxns[j]))
        if τj < τ
            τ = τj
            μ = j
        end
    end

    if t > T return τ end

    μ > 0 ? fire_reaction!(m, reaction(rxns, μ)) : error("No reaction occurred!")

    return τ
end
