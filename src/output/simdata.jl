immutable SimData
    rlz         :: Vector{Realization}
    id2ind      :: Dict{Symbol,Int}
    is_prealloc :: Bool
    metadata    :: Dict{Symbol,UTF8String}

    function SimData(rlz::Vector{Realization}, id2ind::Dict{Symbol,Int}, is_prealloc::Bool)
        f = is_prealloc ? partial_update! : full_update!
        new(rlz, id2ind, is_prealloc, Dict{Symbol,UTF8String}())
    end
end

function SimData(Xt::Vector{Int}, nrlz::Int, id2ind::Dict{Symbol,Int})
    rlz = Realization[ Realization(Xt) for i in 1:nrlz ]

    return SimData(rlz, id2ind, false)
end

function SimData(Xt::Vector{Int}, nrlz::Int, id2ind::Dict{Symbol,Int}, nstates::Int, t::Float64)
    t   = collect(linspace(0.0, t, nstates))
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

function show(io::IO, s::SimData)
    mdata = s.metadata

    print(io, "[algorithm] ", mdata[:algorithm], "\n")
    print(io, " * end time: ",   mdata[:time], "\n")
    print(io, " * no. rlz:  ", length(s.rlz), "\n")
    print(io, " * output:   ", mdata[:datatype])
end

function compile_metadata!(Xt_history::SimData, algorithm::Algorithm, nrlz::Integer, datatype::OutputType)
    mdata = Xt_history.metadata

    alg = string(typeof(algorithm))

    mdata[:algorithm] = utf8(replace(alg, "BioSimulator.", "")) # move this to algorithm constructors
    mdata[:time] = utf8(string(end_time(algorithm)))
    mdata[:nrlz] = utf8(string(nrlz))
    mdata[:datatype] = utf8(isa(datatype, Uniform) ? "partial" : "full")

    tags = get_tags(algorithm)

    for tag in tags
        mdata[tag] = utf8(string(getfield(algorithm, tag)))
    end

    return Xt_history
end
