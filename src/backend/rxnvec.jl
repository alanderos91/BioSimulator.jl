@generated function compute_mass_action(x::Vector{Int}, u::Vector{Int})
    ex = quote
        value = 1
        @inbounds for i in 1:length(u)
            if u[i] > 0
                value = value * x[i]
                @inbounds for k in 2:u[i]
                    value = value * ( x[i] - (k-1) )
                end
            end
        end
        value
    end

    return ex
end

@generated function compute_mass_action(x::Vector{Int}, u::SparseMatrixCSC{Int,Int})
    ex = quote
        value = 1
        @inbounds for i in 1:length(u)
            if u[i] > 0
                value = value * x[i]
                @inbounds for k in 2:u[i]
                    value = value * ( x[i] - (k-1) )
                end
            end
        end
        value
    end

    return ex
end

type PropensityVector{T<:AbstractFloat}
    a  :: Vector{T}
    a0 :: T
end

function PropensityVector{T<:AbstractFloat}(::Type{T}, n)
    return PropensityVector(zeros(T, n), convert(T, Inf))
end

length{T}(pv::PropensityVector{T}) = length(pv.a)
getindex{T}(pv::PropensityVector{T}, i) = getindex(pv.a, i)
setindex!{T}(pv::PropensityVector{T}, val, i) = setindex!(pv.a, val, i)
eachindex{T}(pv::PropensityVector{T}) = eachindex(pv.a)
intensity{T}(pv::PropensityVector{T}) = pv.a0

abstract AbstractReactionSystem

increments(rs::AbstractReactionSystem)   = rs.v
reactants(rs::AbstractReactionSystem)    = rs.u
propensities(rs::AbstractReactionSystem) = rs.a
rates(rs::AbstractReactionSystem)        = rs.k

function compute_propensities!(rs::AbstractReactionSystem, Xt, p)
    a = propensities(rs)
    k = rates(rs)
    a0 = 0.0
    @inbounds for j in eachindex(a)
        key  = k[j]
        kj   = p[key].value
        a[j] = rs(Xt, kj, j)
        a0   = a0 + a[j]
    end
    setfield!(a, :a0, a0)
    return a0
end

immutable SparseReactionSystem <: AbstractReactionSystem
    v :: SparseMatrixCSC{Int,Int}
    u :: SparseMatrixCSC{Int,Int}
    a :: PropensityVector{Float64}
    k :: Vector{Symbol}

    function SparseReactionSystem(v, u, k)
        a = PropensityVector(Float64, length(k))

        new(v, u, a, k)
    end
end

function SparseReactionSystem(dict, id2ind, c, d)
    v = spzeros(c,d)
    u = spzeros(c,d)
    k = Array(Symbol, d)

    j = 1
    for (key, r) in dict
        k[j] = r.rate
        reactants = r.reactants
        products  = r.products

        for (id, i) in id2ind
            x = get(reactants, id, 0)
            y = get(products,  id, 0)
            v[i,j] = y - x
            u[i,j] = x
        end

        j = j + 1
    end

    return SparseReactionSystem(convert(SparseMatrixCSC{Int,Int}, v),
    convert(SparseMatrixCSC{Int,Int}, u),
    k)
end

function fire_reaction!(Xt, rs::SparseReactionSystem, j)
    v    = increments(rs)
    vj   = nonzeros(v)
    idxs = rowvals(v)
    @inbounds for k in nzrange(v, j)
        i = idxs[k]
        Xt[i] = Xt[i] + vj[k]
    end
end

immutable DenseReactionSystem <: AbstractReactionSystem
    v    :: Vector{Vector{Int}}
    u    :: Vector{Vector{Int}}
    a    :: PropensityVector{Float64}
    k    :: Vector{Symbol}

    function DenseReactionSystem(v, u, k)
        inc  = Vector[ v[:, j] for j in 1:size(v, 2) ]
        rct  = Vector[ u[:, j] for j in 1:size(u, 2) ]
        a    = PropensityVector(Float64, length(k))

        new(inc, rct, a, k)
    end
end

function DenseReactionSystem(dict, id2ind, c, d)
    v = Array(Int, c,d)
    u = Array(Int, c,d)
    k = Array(Symbol, d)

    j = 1
    for (key, r) in dict
        k[j] = r.rate
        reactants = r.reactants
        products  = r.products

        for (id, i) in id2ind
            x = get(reactants, id, 0)
            y = get(products,  id, 0)
            v[i,j] = y - x
            u[i,j] = x
        end

        j = j + 1
    end

    return DenseReactionSystem(v, u, k)
end

# propensity of reaction j, assuming rate constant k
function call(rs::DenseReactionSystem, Xt::Vector{Int}, k::Float64, j::Integer)
    u = reactants(rs)
    return k * compute_mass_action(Xt, u[j])
end

# propensity of reaction j, assuming rate constant k
function call(rs::SparseReactionSystem, Xt::Vector{Int}, k::Float64, j::Integer)
    u = reactants(rs)
    return k * compute_mass_action(Xt, u[:,j])
end

function fire_reaction!(Xt, rs::DenseReactionSystem, j)
    v    = increments(rs)
    vj   = v[j]
    @inbounds for k in eachindex(vj)
        Xt[k] = Xt[k] + vj[k]
    end
end

function reaction_system(dict, id2ind)
    d = length(dict)   # number of reactions
    c = length(id2ind) # number of species

    if d < 10
        rs = DenseReactionSystem(dict, id2ind, c, d)
    else
        rs = SparseReactionSystem(dict, id2ind, c, d)
    end

    return rs
end

function mass_action_deriv(Xt, rs, p, i, k)
    key = rates(rs)[i]
    param = p[key]
    mass_action_deriv(Xt, rs, param.value, i, k)
end

# derivative of reaction i with respect to reactant k
function mass_action_deriv(Xt, rs::DenseReactionSystem, rate::Float64, i, k)
    u   = reactants(rs)
    ui  = u[i]
    val = ui[k]

    acc = 0.0
    if val == 0
        acc = 0.0
    elseif val == 1
        acc = helper1(Xt, ui, k)
    elseif val == 2
        acc = 2.0 * Xt[k] - 1
    elseif val > 2
        error("higher order reactions not supported.")
    end

    return rate * acc
end

function mass_action_deriv(Xt, rs::SparseReactionSystem, rate::Float64, i, k)
    u = reactants(rs)
    val = u[k,i]

    acc = 0.0
    if val == 0
        acc = 0.0
    elseif val == 1
        acc = helper2(Xt, u, i, k)
    elseif val == 2
        acc = 2.0 * Xt[k] - 1
    elseif val > 2
        error("higher order reactions not supported.")
    end

    return rate * acc
end

function helper1(Xt::Vector{Int}, u::Vector{Int}, k::Int)
    acc = 1.0
    for i in eachindex(u)
        if i != k
            c = u[i]
            for j in 1:c
                acc = acc * (Xt[i] - (j - 1))
            end
        end
    end
    return acc
end

function helper2(Xt, u, l, k)
    participants = rowvals(u)
    coeff = nonzeros(u)
    a = 1.0
    for i in nzrange(u, l)
        if i != k + length(nzrange(u, l)) - 1
            c = coeff[i]
            ind = participants[i]
            for n in 1:c
                a = a * (Xt[ind] - (n-1))
            end
        end
    end
    return a
end

function fire_reactions!(Xt, rs, events)
    for j in eachindex(events)
        fire_reaction!(Xt, rs, j, events[j])
    end
end

function fire_reaction!(Xt, rs::SparseReactionSystem, j, n)
    v = increments(rs)
    vj = nonzeros(v)
    idxs = rowvals(v)
    for k in nzrange(v, j)
        i = idxs[k]
        Xt[i] = Xt[i] + n * vj[k]
    end
end

function fire_reaction!(Xt, rs::DenseReactionSystem, j, n)
    v = increments(rs)
    vj = v[j]
    for k in eachindex(vj)
        Xt[k] = Xt[k] + n * vj[k]
    end
end
