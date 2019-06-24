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

function Base.show(io :: IO, x :: State)
  println(io, "  type: ", x.ptype)
  print(io, "  nbclass: ", x.class)
end

struct Site{D,T}
  label :: Int
  state :: State
  coord :: MVector{D,T}

  function Site{D,T}(id :: Integer, state :: State, coord :: NTuple{D,T}) where {D,T}
      new{D,T}(id, state, MVector(coord...))
  end
end

Site(id, state, coord::NTuple{N,T}) where {N,T} = Site{N,T}(id, state, coord)

# accessors #

label(x :: Site) = x.label
state(x :: Site) = x.state
coordinates(x :: Site) = x.coord

# for search / sort compatibility
label(x::Int) = x
coordinates(x::NTuple) = MVector(x)
coordinates(x::AbstractVector) = x

get_ptype(x :: Site) = get_ptype(state(x))
get_neighbor_class(x :: Site) = get_neighbor_class(state(x))

# state #

transform!(x :: Site, l) = transform!(state(x), l)
change_neighbor_class!(x :: Site, k) = change_neighbor_class!(state(x), k)

# iteration #
Base.iterate(x::Site) = iterate(coordinates(x))
Base.iterate(x::Site, state) = iterate(coordinates(x), state)

# indexing #

Base.getindex(x :: Site, i) = getindex(coordinates(x), i)
Base.setindex!(x :: Site, val, i) = setindex!(coordinates(x), val, i)
Base.eachindex(x :: Site) = eachindex(coordinates(x))

Base.isless(x::Site, y::Site) = isless(label(x), label(y))

# IO #

Base.summary(io::IO, x::Site) = print(io, "Site")

function Base.show(io::IO, x::Site)
  # summary(io, x)
  print(io, coordinates(x))
end

function Base.show(io::IO, m::MIME"text/plain", x::Site)
  lk = string(get_ptype(x), "|", get_neighbor_class(x))
  print(io, lk)
  print(io, " @ ")
  print(io, coordinates(x))
end

##### recipes #####

@recipe function f(::Type{Site{1,T}}, x::Site{1,T}) where T
  [(0, x[1])]
end

@recipe function f(::Type{Site{2,T}}, x::Site{2,T}) where T
  [(x[1], x[2])]
end

@recipe function f(::Type{Site{3,T}}, x::Site{3,T}) where T
  [(x[1], x[2], x[3])]
end

function build_composition!(comp, lattice, x::Site{D}) where D
  fill!(comp, zero(eltype(comp)))
  comp[1] = capacity(topology(lattice), D)

  for id in neighborhood(lattice, x)
    y = get_site(lattice, id)
    l = get_ptype(y)

    comp[1] -= 1
    comp[l] += 1
  end

  return comp
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
