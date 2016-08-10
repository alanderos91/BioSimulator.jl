immutable History
    history  :: Vector{Int}
    target   :: Int

    function History(h::Vector{Int}, target::Int)
        new(h, target)
    end
end

# no prealloc
function History(target::Int)
    return History(Int[], target)
end

# yes prealloc
function History(target::Int, nstates::Int)
    return History(zeros(Int, nstates), target)
end

function partial_update!(h::History, Xt::Vector{Int}, n::Int)
    hist = h.history
    i    = h.target

    hist[n] = Xt[i]

    return nothing
end

function full_update!(h::History, Xt::Vector{Int}, n::Int)
    hist = h.history
    i    = h.target

    push!(hist, Xt[i])

    return nothing
end

getindex(h::History, i) = h.history[i]

function Base.show(io::IO, h::History)
    show(io, h.history)
end
