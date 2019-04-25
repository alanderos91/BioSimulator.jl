mutable struct State
  ptype :: Int
  class :: Int
end

# outer constructors

State() = State(1, 1)
State(ptype) = State(ptype, 1)

get_ptype(x :: State) = x.ptype
get_neighbor_class(x :: State) = x.class

transform!(x :: State, l) = (x.ptype = l)
change_neighbor_class!(x :: State, k) = (x.class = k)

show(io :: IO, x :: State) = print(io, "(l, k) = ", x.ptype, ", ", x.class)

struct Site{D,T,M}
  label :: Int
  state :: State
  coord :: MVector{D,T}
  neigh :: M

  function Site{D,T,M}(id :: Integer, state :: State, coord :: NTuple{D,T}) where {D,T,M}
      new{D,T,M}(id, state, MVector(coord...), M())
  end
end

Site(id, state, coord::NTuple{1,T}) where T = Site{1,T,VonNeumann1D{Int}}(id, state, coord)
Site(id, state, coord::NTuple{2,T}) where T = Site{2,T,VonNeumann2D{Int}}(id, state, coord)
Site(id, state, coord::NTuple{3,T}) where T = Site{3,T,VonNeumann3D{Int}}(id, state, coord)

# Site(id, state, coord, Hexagonal2D)
Site(id, state, coord::NTuple{1,T}, ::Type{M}) where {T, M} = Site{1,T,M}(id, state, coord)
Site(id, state, coord::NTuple{2,T}, ::Type{M}) where {T, M} = Site{2,T,M}(id, state, coord)
Site(id, state, coord::NTuple{3,T}, ::Type{M}) where {T, M} = Site{3,T,M}(id, state, coord)

# accessors #

label(x :: Site) = x.label
state(x :: Site) = x.state
coordinates(x :: Site) = x.coord

get_ptype(x :: Site) = get_ptype(state(x))
get_neighbor_class(x :: Site) = get_neighbor_class(state(x))

neighborhood(x :: Site) = x.neigh
neighborcount(x :: Site) = neighborcount(neighborhood(x))

# state #

transform!(x :: Site, l) = transform!(state(x), l)
change_neighbor_class!(x :: Site, k) = change_neighbor_class!(state(x), k)

# iteration #
iterate(x::Site) = iterate(coordinates(x))
iterate(x::Site, state) = iterate(coordinates(x), state)

# indexing #

getindex(x :: Site, i) = getindex(coordinates(x), i)
setindex!(x :: Site, val, i) = setindex!(coordinates(x), val, i)
eachindex(x :: Site) = eachindex(coordinates(x))

# v0.6 specific
import Base: isless
Base.isless(x::Site, y::Site) = isless(label(x), label(y))

# operations #

# # candidate for inlining; generated function?
# function swap_sites!(x :: Site{D}, y :: Site{D}) where D
#     for i in eachindex(x)
#         x[i], y[i] = y[i], x[i]
#     end
#     return nothing
# end

eachdir(x :: Site{D,T,M}) where {D,T,M} = eachdir(M, x)

# IO #

show(io :: IO, x :: Site) = print(io, label(x), " / ", state(x), " / ", coordinates(x))

# neighborhood #

add_neighbor!(x :: Site, y :: Site) = push!(neighborhood(x), y)
rmv_neighbor!(x :: Site, y :: Site) = delete!(neighborhood(x), y)

function build_composition!(composition, x :: Site{D,T,M}) where {D,T,M}
  fill!(composition, zero(eltype(composition)))

  composition[1] = 2 * D

  for y in neighborhood(x)
      l = get_ptype(y)
      composition[1] -= 1
      composition[l] += 1
  end

  return nothing
end

const XD1 = (1,)

const XD2 = (1, 2)
const YD2 = (2, 1)

const XD3 = (1, 2, 3)
const YD3 = (1, 3, 2)
const ZD3 = (2, 3, 1)

function lex_order(x, y, idx)
  for i in idx
    a, b = x[i], y[i]
    if !isequal(a, b)
      return isless(a, b)
    end
  end
  return false
end

lex_order_x1(x, y) = lex_order(x, y, XD1)

lex_order_x2(x, y) = lex_order(x, y, XD2)
lex_order_y2(x, y) = lex_order(x, y, YD2)

lex_order_x3(x, y) = lex_order(x, y, XD3)
lex_order_y3(x, y) = lex_order(x, y, YD3)
lex_order_z3(x, y) = lex_order(x, y, ZD3)

function lex_sort!(v; alg = Base.Sort.DEFAULT_STABLE, lt = lex_order_x, by = identity, rev::Bool = false, order :: Base.Ordering = Base.Order.Forward)
  sort!(v, alg = alg, lt = lt, by = by, rev = rev, order = order)
end

@generated function distance(
  x :: Site{D,T,M},
  y :: Site{D,T,M}) where {D,T,M}
  ex = :(0)

  for i = 1:D
      ex = :($ex + abs(x[$i] - y[$i]))
  end

  return ex
end

function build_local_neighborhoods!(sites::Vector{Site{1,T,M}}) where {T,M}
  lex_sort!(sites, lt = lex_order_x1, by = coordinates)
  sweep_neighbors!(sites)

  return nothing
end

function build_local_neighborhoods!(sites::Vector{Site{2,T,M}}) where {T,M}
  lex_sort!(sites, lt = lex_order_x2, by = coordinates)
  sweep_neighbors!(sites)
  
  lex_sort!(sites, lt = lex_order_y2, by = coordinates)
  sweep_neighbors!(sites)

  return nothing
end

function build_local_neighborhoods!(sites::Vector{Site{3,T,M}}) where {T,M}
  lex_sort!(sites, lt = lex_order_x3, by = coordinates)
  sweep_neighbors!(sites)
  
  lex_sort!(sites, lt = lex_order_y3, by = coordinates)
  sweep_neighbors!(sites)
  
  lex_sort!(sites, lt = lex_order_z3, by = coordinates)
  sweep_neighbors!(sites)

  return nothing
end

function sweep_neighbors!(sites)
  for i = 2:length(sites)
    x = sites[i - 1]
    y = sites[i]
    d = distance(x, y)
    if d == 1
      add_neighbor!(x, y)
      add_neighbor!(y, x)
    end
  end
  return sites
end
