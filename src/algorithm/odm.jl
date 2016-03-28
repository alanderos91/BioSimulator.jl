type ODM <: Algorithm
    # parameters
    T::Float64
    init_steps::Int
    init_iters::Int
    coupled::Symbol
    callback::Function

    # state variables
    t::Float64
    steps::Int
    g::DiGraph

    # statistics
    avg_nsteps::Float64
    avg_step_size::Float64

    # metadata tags
    tags::Vector{Symbol}

    function ODM(T, init_steps, init_iters, coupled, callback)
        new(T, init_steps, init_iters, coupled, callback, 0.0, 0, DiGraph(), 0.0, 0.0, DEFAULT_TAGS)
    end
end

function odm(T; init_steps=100, init_iters=1, coupled=:tight, na...)
    return ODM(T, init_steps, init_iters, coupled, update_all!)
end

function initialize!(x::ODM, m::Model)
    if x.coupled == :loose
        g = init_dep_graph()
        setfield!(x, :g, g)
        setfield!(x, :callback, update_dependents!)
    end
    # Presimulate to sort reactions according to multiscale property.
    # This will modify spcs and rxns
    presimulate!(x, m)
    return;
end

function reset!(x::ODM, m::Model)
    setfield!(x, :t, 0.0)
    setfield!(x, :steps, 0)
    compute_propensities!(m)

    return;
end

function call(x::ODM, m)
    τ = odm_update!(m ,x)

    # update algorithm variables
    setfield!(x, :t,     time(x) + τ)
    setfield!(x, :steps, steps(x) + 1)

    compute_statistics!(x, τ)
end

function odm_update!(m::Model, x)
    t = x.t
    T = x.T

    rxn = reactions(m)
    a0 = intensity(rxn)

    τ = rand(Exponential(1 / a0))

    if t + τ > T return τ end

    μ = select_reaction(m, a0)
    fire_reaction!(m, reaction(rxn, μ))

    update_propensities = getfield(x, :callback)
    update_propensities(x, m, μ)

    return τ
end

function update_dependents!(x::ODM, m::Model, μ)
    g = getfield(x, :g)
    Xt = species(m)
    rxns = reactions(m)
    p = parameters(m)
    a0 = intensity(rxns)
    dependents = neighbors(g, μ)

    @inbounds for α in dependents
        a0 = a0 - rxns[α]
        rxns[α] = intensity(reaction(rxns, α), Xt, p)
        a0 = a0 + rxns[α]
    end
    setfield!(rxn, :intensity, a0)
    return;
end

function update_all!(x::ODM, m::Model, μ)
    compute_propensities!(m)
    return;
end

function presimulate!(x, m)
    itr = x.init_iters
    n   = x.init_steps

    rxns = reactions(m)
    events = zeros(Float64, length(rxns))

    for i = 1:itr
        reset!(m)
        for k = 1:n
            a0 = compute_propensities!(m)
            τ = rand(Exponential(1 / a0))
            jump = a0 * rand()
            μ = select_reaction(m, a0)
            fire_reaction!(m, reaction(rxns, μ))
            events[μ] = events[μ] + 1
        end
    end

    for i in eachindex(rxns)
        @inbounds rxns[i] = events[i] / itr
    end

    sort!(rxns)
end
