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

function call(x::SSA, m::Model)
    compute_propensities!(m)
    τ = ssa_update!(m, x)

    # update algorithm variables
    setfield!(x, :t,     time(x) + τ)
    setfield!(x, :steps, steps(x) + 1)

    compute_statistics!(x, τ)
end

function ssa_update!(m::Model, x)
    t = x.t
    T = x.T
    rxn = reactions(m)
    a0 = intensity(rxn)

    τ = rand(Exponential(1 / a0))

    if t + τ > T return τ end

    μ = select_reaction(m, a0)
    fire_reaction!(m, reaction(rxn, μ))

    return τ
end

function select_reaction(m::Model, a0)
    jump = a0 * rand()
    rxn = reactions(m)
    s = 0.0

    for i in eachindex(rxn)
        s = s + rxn[i]
        if s >= jump
            return i
        end
    end
    error("no reaction occurred!")
end
