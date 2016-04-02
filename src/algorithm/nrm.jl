import Base.Collections: PriorityQueue, peek
import Base.Order: ForwardOrdering

type NRM <: Algorithm
    # parameters
    T::Float64

    # state variables
    t::Float64
    g::DiGraph
    pq::PriorityQueue{Int,Float64,ForwardOrdering}
    propensities::Vector{Float64}
    steps::Int

    # statistics
    avg_nsteps::Float64
    avg_step_size::Float64

    # metadata tags
    tags::Vector{Symbol}

    function NRM(T)
        new(T, 0.0, DiGraph(), PriorityQueue(Int,Float64), Float64[], 0, 0.0, 0.0, DEFAULT_TAGS)
    end
end

function nrm(T; na...)
    return NRM(T)
end

function initialize!(x::NRM, m)
    rxns = reactions(m)
    d    = length(rxns)

    g = init_dep_graph(rxns)
    setfield!(x, :g,  g)

    compute_propensities!(m)

    a = propensities(rxns)
    setfield!(x, :propensities, deepcopy(a))

    pq = init_pq(x)
    setfield!(x, :pq, pq)
    return;
end

function reset!(x::NRM, m::Model)
    setfield!(x, :t, 0.0)
    setfield!(x, :steps, 0)
    compute_propensities!(m)
    init_pq!(x)

    return;
end

function step!(x::NRM, m::Model)
    τ = nrm_update!(m ,x)

    # update algorithm variables
    setfield!(x, :t,     τ)
    setfield!(x, :steps, steps(x) + 1)

    update_dependents!(x, m)
    compute_statistics!(x, τ)
end

function nrm_update!(m::Model, x)
    T = end_time(x)
    rxns = reactions(m)
    pq = getfield(x, :pq)
    μ, τ = peek(pq)

    if τ > T return τ end

    μ > 0 ? fire_reaction!(m, reaction(rxns, μ)) : error("No reaction ocurred!")
    return τ
end

function init_dep_graph(rxns::ReactionVector)
    d = length(rxns)
    g = DiGraph(d)

    for j in eachindex(rxns)
        r = reaction(rxns, j)
        pre = r.pre
        post = r.post

        # Search for reactions r_i with reactant or product of r_j as reactant
        # how to handle immigration reaction?
        for k in eachindex(pre)
            if pre[k] != 0 || post[k] != 0
                for i in eachindex(rxns)
                    if i == j continue end
                    if reaction(rxns, i).pre[k] != 0 add_edge!(g, j, i) end
                end
            end
        end
    end
    return g
end

function update_dependents!(x::NRM, m::Model)
    pq   = getfield(x, :pq)
    g    = getfield(x, :g)
    μ, τ = peek(pq)

    t    = time(x)
    Xt   = species(m)
    rxns = reactions(m)
    p    = parameters(m)

    # Compute propensities and firing times for dependent reactions
    dependents = neighbors(g, μ)
    for α in dependents
        old_propensity = rxns[α]
        rxns[α] = intensity(reaction(rxns, α), Xt, p)
        pq[α] = t + (old_propensity / rxns[α]) * (pq[α] - t)
    end
    rxns[μ] = intensity(reaction(rxns, μ), Xt, p)
    pq[μ]   = τ + rand(Exponential(1 / rxns[μ]))
    return g
end

function init_pq(x)
    a  = getfield(x, :propensities)
    d  = length(a)
    pq = PriorityQueue(collect(1:d), zeros(Float64, d))
    return pq
end

function init_pq!(x)
    a  = getfield(x, :propensities)
    pq = getfield(x, :pq)
    for j in eachindex(a)
        pq[j] = rand(Exponential(1 / a[j]))
    end
end
