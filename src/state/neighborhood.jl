abstract type Neighborhood end

# accessors #

neighbors(x :: Neighborhood) = x.neighbors
capacity(x :: Neighborhood) = x.capacity

neighborcount(x :: Neighborhood) = length(neighbors(x))

# Base overloads #

function push!(x :: Neighborhood, y :: Site)
    neighborcount(x) ≥ capacity(x) && throw(ErrorException(string("neighborhood at capacity", x, "\n", "neighbor = ", y)))
    push!(neighbors(x), y)
end

function delete!(x :: Neighborhood, y :: Site)
    y ∉ x && throw(ErrorException("attempted to remove site not in neighborhood."))
    ix = findfirst(x, y)
    deleteat!(neighbors(x), ix)
end

# iteration #
iterate(x::Neighborhood) = iterate(neighbors(x))
iterate(x::Neighborhood, state) = iterate(neighbors(x), state)

eachindex(x :: Neighborhood) = eachindex(neighbors(x))
length(x :: Neighborhood) = length(neighbors(x))

getindex(x :: Neighborhood, inds...) = getindex(neighbors(x), inds...)
setindex!(x :: Neighborhood, vals, inds...) = setindex!(neighbors(x), vals, inds...)

in(x :: Site, m :: Neighborhood) = in(x, neighbors(m))
empty!(x :: Neighborhood) = empty!(neighbors(x))

# sampling #

function sample(x :: Neighborhood, l, nl)
    c = nl * rand()
    j = 1
    s = get_ptype(x[j]) == l ? 1 : 0

    while s < c && j < length(x)
        isequal(get_ptype(x[j += 1]), l) && (s += 1)
    end

    return x[j]
end
