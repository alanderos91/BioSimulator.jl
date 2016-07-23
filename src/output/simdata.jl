immutable SimData
    rlz      :: Vector{Realization}
    id2ind   :: Dict{Symbol,Int}
    is_prealloc :: Bool

    function SimData(rlz::Vector{Realization}, id2ind::Dict{Symbol,Int}, is_prealloc::Bool)
        f = is_prealloc ? partial_update! : full_update!
        new(rlz, id2ind, is_prealloc)
    end
end

function SimData(Xt::Vector{Int}, nrlz::Int, id2ind::Dict{Symbol,Int})
    rlz = Realization[ Realization(Xt) for i in 1:nrlz ]

    return SimData(rlz, id2ind, false)
end

function SimData(Xt::Vector{Int}, nrlz::Int, id2ind::Dict{Symbol,Int}, nstates::Int)
    t   = zeros(Float64, nstates)
    rlz = Realization[ Realization(Xt, t) for i in 1:nrlz ]

    return SimData(rlz, id2ind, true)
end

function update!(s::SimData, i::Int, t::Float64, Xt::Vector{Int})
    r = s.rlz[i]
    is_prealloc = s.is_prealloc

    if is_prealloc
        partial_update!(r, t, Xt)
    else
        full_update!(r, t, Xt)
    end

    return nothing
end

function show(s::SimData)
    println("n. realization", length(s.rlz))
end
