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
        a[j] = rs(Xt, kj, j)::Float64
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
    klaw :: Vector{Function}

    function SparseReactionSystem(v, u, k)
        a = PropensityVector(Float64, length(k))
        klaw = make_propensities(u)
        new(v, u, a, k, klaw)
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
    klaw :: Vector{Function}

    function DenseReactionSystem(v, u, k)
        inc  = Vector[ v[:, j] for j in 1:size(v, 2) ]
        rct  = Vector[ u[:, j] for j in 1:size(u, 2) ]
        a    = PropensityVector(Float64, length(k))
        klaw = make_propensities(u)

        new(inc, rct, a, k, klaw)
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

function make_propensity(u, id)
    name = symbol("_a", id)
    args = Expr[]

    for i in 1:length(u)
        if u[i] > 0
            push!(args, :( x[$i] ))
            for k in 2:u[i]
                push!(args, :( (x[$i] - $(k-1)) ))
            end
        end
    end

    ex   = quote end
    code = Expr(:(=), :( $name(x::Vector{Int}) ))

    if length(args) > 1
        ex   = Expr(:call, :*, args...)
        push!(code.args, ex)
    elseif length(args) == 1
        ex   = args[1]
        push!(code.args, ex)
    else
        push!(code.args, 1)
    end

    return eval(code)
end

function make_propensities(u)
    klaw = Function[]

    for i in 1:size(u,2)
        f = make_propensity(u[:,i], i)
        push!(klaw, f)
    end

    return klaw
end

# propensity of reaction j, assuming rate constant k
function call(rs::AbstractReactionSystem, Xt::Vector{Int}, k::Float64, j::Integer)
    a = rs.klaw[j]
    return k * a(Xt)::Int
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
