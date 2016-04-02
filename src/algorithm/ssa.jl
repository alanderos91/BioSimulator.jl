type SSA <: Algorithm
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

    function SSA(T)
        new(T, 0.0, 0, 0.0, 0.0, DEFAULT_TAGS)
    end
end

function ssa(T; na...)
    return SSA(T)
end

initialize!(x::SSA, m::Model) = return;

function reset!(x::SSA, m::Model)
    setfield!(x, :t, 0.0)
    setfield!(x, :steps, 0)
    return;
end

function step!(x::SSA, m::Model)
    a0 = compute_propensities!(m)
    τ = rand(Exponential(1 / a0))

    # update algorithm variables
    setfield!(x, :t,     time(x) + τ)
    setfield!(x, :steps, steps(x) + 1)
    compute_statistics!(x, τ)

    if !done(x, m)
        r = reactions(m)
        μ = select_reaction(r, a0)

        if μ > 0
            fire_reaction!(m, reaction(r, μ))
        else
            error("no reaction occurred!")
        end
    end
end


function select_reaction(r::ReactionVector, a0)
    jump = a0 * rand()
    a = 0.0

    for j in eachindex(r)
        a = a + r[j]
        if a >= jump return j end
    end
    return 0
end
