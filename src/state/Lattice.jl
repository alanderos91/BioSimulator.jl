struct Lattice{D,T,M,U}
  site::Vector{Site{D,T}}
  neighbors::Vector{Vector{Int}}
  types::Vector{U}
end

function Lattice(coord::Matrix, types::Vector; nbhood = VonNeumann())
  dimension = size(coord, 1)
  unique_types = unique(types)
  number_types = length(unique_types)
  number_neighbors = capacity(nbhood, dimension)

  labels = Dict(unique_types[i] => i for i in eachindex(unique_types))
  site = [Site(i, State(labels[types[i]]), tuple(coord[:,i]...)) for i in eachindex(types)]

  neighbors = [sizehint!(Int[], number_neighbors) for i in eachindex(site)]

  mapping = [(label => l + 1) for (label, l) in labels]
  T = eltype(coord)
  U = eltype(mapping)

  return Lattice{dimension,T,typeof(nbhood),U}(site, neighbors, mapping)
end

##### accessors

dimension(::Lattice{D}) where D <: Number = D
topology(::Lattice{D,T,M}) where {D,T,M} = M
number_types(x::Lattice) = length(x.types)
number_sites(x::Lattice) = length(x.site)

##### IO
function Base.summary(io::IO, ::Lattice{D,T,M}) where {D,T,M}
  print(io, "$(D)-D Lattice with $(M) neighborhoods")
end

function Base.show(io::IO, lattice::Lattice)
  Base.summary(io, lattice)
  print(io, "\n  species: ")
  show(io, lattice.types)
end

function Base.show(io::IO, m::MIME"text/plain", lattice::Lattice)
  Base.summary(io, lattice)
  print(io, "\n", "species: ")
  show(io, m, lattice.types)
end

##### other Base overloads
Base.copy(lattice::Lattice) = deepcopy(lattice)

##### query API

function istracked(lattice::Lattice, coord)
  # ?
end

# get the site living at the given coordinates
function get_site(lattice::Lattice, coord)
  idx = searchsorted(lattice.site, coord, by = coordinates)

  # if idx is empty, site was not found

  # if idx has length > 1, we have a duplicate

  # if idx has length = 1, return the object
end

# get the site with the given id
function get_site(lattice::Lattice, id::Int)
  idx = searchsorted(lattice.site)

  # if idx is empty, site was not found

  # if idx has length > 1, we have a duplicate

  # if idx has length = 1, return the object
end

function neighborhood(lattice::Lattice, site)
  lattice.nbhood[label(site)]
end

##### update API

function spawn_new_site(lattice::Lattice, coord)
  id = number_sites(lattice) + 1

  return Site(label, State(), coord)
end

function add_site!(lattice::Lattice{D,T,M}, new_site) where {D,T,M}
  # add neighborhood for site
  push!(lattice.neighbors, sizehint!(Int[], capacity(M(), D)))

  # add site to list
  push!(lattice.site, new_site)
end

##### recipes

@recipe function f(lattice::Lattice{D,T,VonNeumann}) where {D,T}
  site = lattice.site
  types = lattice.types

  seriestype := :scatter
  aspect_ratio --> 1
  markerstrokewidth --> 0
  markershape --> :square
  # color --> map(get_ptype, site)
  # label --> map(first, lattice.types)

  for t in types
    myseries = filter(x -> get_ptype(x) == last(t) - 1, site)

    @series begin
      label --> first(t)
      [tuple(s...) for s in myseries]
    end
  end
end

@recipe function f(lattice::Lattice{2,T,Hexagonal}) where {T}
  site = lattice.site
  types = lattice.types

  seriestype := :scatter
  aspect_ratio --> 1
  markerstrokewidth --> 0
  markershape --> :hexagon
  # color --> map(get_ptype, site)
  # label --> map(first, lattice.types)

  A = cos(pi / 3)
  B = sin(pi / 3)

  for t in types
    s = filter(x -> get_ptype(x) == last(t) - 1, site)

    @series begin
      label --> first(t)
      new_coords = Vector{NTuple{2,Float64}}(undef, length(s))

      for i in eachindex(s)
        p = s[i]

        x = p[1] + A * p[2]
        y = B * p[2]

        new_coords[i] = (x, y)
      end

      new_coords
    end
  end
end
