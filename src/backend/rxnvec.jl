@generated function compute_mass_action(x::Vector{Int}, u::Vector{Vector{Int}}, j::Int)
    ex = quote
        r = u[j]
        value = 1

        @inbounds for i in 1:length(r)
            if r[i] > 0
                value = value * x[i]
                @inbounds for k in 2:r[i]
                    value = value * ( x[i] - (k-1) )
                end
            end
        end
        value
    end

    return ex
end

@generated function compute_mass_action(x::Vector{Int}, u::SparseMatrixCSC{Int,Int}, j::Int)
    ex = quote
        rv = rowvals(u)
        nz = nonzeros(u)
        value = 1

        @inbounds for i in nzrange(u, j)
            if nz[i] > 0
                value = value * x[rv[i]]
                @inbounds for k in 2:nz[i]
                    value = value * ( x[rv[i]] - (k-1) )
                end
            end
        end
        value
    end

    return ex
end

@generated function compute_mass_action_deriv(Xt::Vector{Int}, u::Vector{Vector{Int}}, j::Int, k::Int)
    ex = quote
        r = u[j]
        value = 0

        if r[k] > 0
            value = 1
            @inbounds for i in eachindex(r)
                if i != k
                    if r[i] > 0
                        value = value * Xt[i]
                        @inbounds for n in 2:r[i]
                            value = value * (Xt[i] - (n - 1))
                        end
                    end
                elseif r[i] > 1
                    tmp1 = 0
                    @inbounds for m in 1:r[i]
                        tmp2 = 1
                        @inbounds for n in 1:r[i]
                            if m != n; tmp2 = tmp2 * (Xt[i] - (n - 1)); end
                        end
                        tmp1 = tmp1 + tmp2
                    end
                    value = value * tmp1
                end
            end
        end

        value
    end
    return ex
end

@generated function compute_mass_action_deriv(Xt::Vector{Int}, u::SparseMatrixCSC{Int,Int}, j::Int, k::Int)
    ex = quote
        rv = rowvals(u)
        nz = nonzeros(u)
        value = 0

        if u[k,j] > 0
            value = 1
            @inbounds for i in nzrange(u, j)
                if rv[i] != k
                    if nz[i] > 0
                        value = value * Xt[rv[i]]
                        @inbounds for n in 2:nz[i]
                            value = value * (Xt[rv[i]] - (n - 1))
                        end
                    end
                elseif nz[i] > 1
                    tmp1 = 0
                    @inbounds for m in 1:nz[i]
                        tmp2 = 1
                        @inbounds for n in 1:nz[i]
                            if m != n; tmp2 = tmp2 * (Xt[rv[i]] - (n - 1)); end
                        end
                        tmp1 = tmp1 + tmp2
                    end
                    value = value * tmp1
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
function call(rs::AbstractReactionSystem, Xt::Vector{Int}, k::Float64, j::Integer)
    u = reactants(rs)
    return k * compute_mass_action(Xt, u, j)
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

# unpack
function mass_action_deriv(Xt, rs, p, i, k)
    key = rates(rs)[i]
    param = p[key]
    mass_action_deriv(Xt, rs, param.value, i, k)
end

# derivative of reaction j with respect to reactant k
function mass_action_deriv(Xt, rs::AbstractReactionSystem, rate::Float64, j, k)
    return rate * compute_mass_action_deriv(Xt, reactants(rs), j, k)
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
