import Base.Collections: PriorityQueue, peek
import Base.Order: ForwardOrdering

"""
```
NRM(T)
```

Gibson and Bruck's Next Reaction Method, equivalent to the original `SSA`.

### Arguments
- `T`: The simulation end time.
"""
type NRM <: ExactMethod
    # parameters
    T :: Float64

    # state variables
    t            :: Float64
    g            :: DiGraph
    pq           :: PriorityQueue{Int,Float64,ForwardOrdering}
    steps        :: Int

    # statistics
    avg_nsteps    :: Float64
    avg_step_size :: Float64

    # metadata tags
    tags :: Vector{Symbol}

    function NRM(T)
        new(T, 0.0, DiGraph(), PriorityQueue(Int,Float64), 0, 0.0, 0.0, DEFAULT_EXACT)
    end
end

function nrm(T; na...)
    return NRM(T)
end

function initialize!(x::NRM, m::Model)
    Xt = species(m)
    rs = reactions(m)
    p  = parameters(m)
    a  = propensities(rs)

    # initialize dependency graph
    g = init_dep_graph(rs)

    setfield!(x, :g,  g)

    # initialize the priority queue
    pq = init_pq(a)
    setfield!(x, :pq, pq)
    return;
end

function reset!(x::NRM, m::Model)
    pq = getfield(x, :pq)

    Xt = species(m)
    rs = reactions(m)
    p  = parameters(m)
    a  = propensities(rs)

    setfield!(x, :t, 0.0)
    setfield!(x, :steps, 0)
    compute_propensities!(rs, Xt, p)
    init_pq!(pq, a)

    return;
end

function step!(x::NRM, Xt, rs, p)
    pq   = getfield(x, :pq)
    μ, τ = peek(pq)
    a0   = intensity(propensities(rs))

    # update algorithm variables
    t = time(x)
    setfield!(x, :t, τ)
    setfield!(x, :steps, steps(x) + 1)
    compute_statistics!(x, τ - t)

    if time(x) < end_time(x) && a0 > 0
        fire_reaction!(Xt, rs, μ)
        update_dependents!(x, Xt, rs, p)
    end

    return;
end

function init_dep_graph(rs::AbstractReactionSystem)
    d = length(rs.a)
    g = DiGraph(d)
    u = reactants(rs)
    v = increments(rs)

    g = populate!(g, u, v, d)
    return g
end

function populate!(g, u::Vector{Vector{Int}}, v::Vector{Vector{Int}}, d)

    for j in 1:d
        u1 = u[j]
        v1 = v[j]
        for k in eachindex(u1)
            if u1[k] != 0 || v1[k] + u1[k] != 0
                for i in 1:d
                    if i == j continue end
                    u2 = u[i]
                    if u2[k] != 0 add_edge!(g, j, i) end
                end
            end
        end
    end
    return g
end

function populate!(g, u::SparseMatrixCSC{Int,Int}, v::SparseMatrixCSC{Int,Int}, d)
    urow = rowvals(u)
    uval = nonzeros(v)
    vrow = rowvals(v)
    vval = nonzeros(v)
    c, d = size(u)

    for j in 1:d
        for k in 1:c
            if u[k,j] != 0 || v[k,j] + u[k,j] != 0
                for i in 1:d
                    if i == j continue end
                    if u[k,i] != 0 add_edge!(g, j, i) end
                end
            end
        end
    end
    return g
end

function update_dependents!(x::NRM, Xt, rs, p)
    pq   = getfield(x, :pq)
    g    = getfield(x, :g)
    μ, τ = peek(pq)
    t    = time(x)
    a    = propensities(rs)
    k    = rates(rs)

    # Compute propensities and firing times for dependent reactions
    dependents = neighbors(g, μ)
    for α in dependents
        key   = k[α]
        kα    = p[key].value
        old_a = a[α]
        a[α]  = rs(Xt, kα, α)
        a.a0  = a.a0 - old_a + a[α]
        if pq[α] != Inf
            pq[α] = t + (old_a / a[α]) * (pq[α] - t)
        else
            pq[α] = τ + rand(Exponential(1 / a[α]))
        end
    end
    key   = k[μ]
    kμ    = p[key].value
    old_a = a[μ]
    a[μ]  = rs(Xt, kμ, μ)
    a.a0  = a.a0 - old_a + a[μ]
    pq[μ] = τ + rand(Exponential(1 / a[μ]))
    return g
end

function init_pq(a)
    d  = length(a)
    pq = PriorityQueue(collect(1:d), zeros(Float64, d))
    init_pq!(pq, a)
end

function init_pq!(pq, a)
    for j in eachindex(a)
        pq[j] = rand(Exponential(1 / a[j]))
    end
    return pq
end
