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

function step!(x::FRM, m::Model)
    compute_propensities!(m)

    r = reactions(m)

    τ, μ = select_reaction(r)

    # update algorithm variables
    setfield!(x, :t,     time(x) + τ)
    setfield!(x, :steps, steps(x) + 1)

    if !done(x, m)
        if μ > 0
            fire_reaction!(m, reaction(r, μ))
        else
            error("no reaction ocurred!")
        end
    end

    compute_statistics!(x, τ)
end

function select_reaction(r::ReactionVector)
    τ = Inf
    μ = 0
    for j in eachindex(r)
        τj = rand(Exponential(1 / r[j]))
        if τj < τ
            τ = τj
            μ = j
        end
    end
    return τ, μ
end
