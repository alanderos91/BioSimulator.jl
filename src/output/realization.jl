type TickTock
    val :: Int

    function TickTock()
        return new(1)
    end
end

function next!(tt::TickTock)
    tt.val += 1
    return tt.val
end

immutable Realization
    t    :: Vector{Float64}
    hist :: Vector{History}
    pos  :: TickTock

    function Realization(t::Vector{Float64}, hist::Vector{History})
        new(t, hist, TickTock())
    end
end

# no prealloc
function Realization(Xt::Vector{Int})
    t    = Float64[]
    hist = History[ History(i) for i in 1:length(Xt) ]
    return Realization(t, hist)
end

# yes prealloc
function Realization(Xt::Vector{Int}, t::Vector{Float64})
    hist = History[ History(i, length(t)) for i in 1:length(Xt) ]
    return Realization(t, hist)
end

function partial_update!(r::Realization, t::Float64, Xt::Vector{Int})
    len = length(r.t)
    pos = r.pos.val
    ttt = r.t

    while (pos <= len) && (t >= ttt[pos])
        @inbounds for h in r.hist
            partial_update!(h, Xt, pos)
        end
        pos = next!(r.pos)
    end

    return nothing
end

function full_update!(r::Realization, t::Float64, Xt::Vector{Int})
    push!(r.t, t)

    @inbounds for h in r.hist
        full_update!(h, Xt, 0)
    end

    return nothing
end

# entire history k
getindex(r::Realization, k) = r.hist[k]

# slice i of history k
getindex(r::Realization, k, i) = r.hist[k][i]
