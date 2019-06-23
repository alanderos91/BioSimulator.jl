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

##### query API

function istracked(lattice::Lattice, coord)
  # ?
end

function get_site(lattice::Lattice, coord)
  idx = searchsorted(lattiec.site, coord, by = coordinates)

  # if idx is empty, site was not found

  # if idx has length > 1, we have a duplicate

  # if idx has length = 1, return the object
end

function neighborhood(lattice::Lattice, site)
  lattice.nbhood[label(site)]
end

##### update API

function add_site!(lattice::Lattice{D,T,M}, new_site) where {D,T,M}
  # add neighborhood for site
  push!(lattice.neighbors, sizehint!(Int[], capacity(M(), D)))

  # add site to list
  push!(lattice.site, new_site)
end

##### recipes

@recipe function f(::Type{Lattice{D,T,VonNeumann}}, lattice::Lattice{D,T,VonNeumann}) where {D,T}
  site = lattice.site

  seriestype := :scatter
  aspect_ratio --> 1
  markerstrokewidth --> 0
  markershape --> :square
  color --> map(get_ptype, site)

  [tuple(s...) for s in site]
end

@recipe function f(::Type{Lattice{2,T,Hexagonal}}, lattice::Lattice{2,T,Hexagonal}) where {T}
  site = lattice.site

  seriestype := :scatter
  aspect_ratio --> 1
  markerstrokewidth --> 0
  markershape --> :hexagon
  color --> map(get_ptype, site)

  new_coords = Vector{NTuple{2,Float64}}(undef, length(site))

  A = cos(pi / 3)
  B = sin(pi / 3)

  for i in eachindex(site)
    s = site[i]

    x = s[1] + A * s[2]
    y = B * s[2]

    new_coords[i] = (x, y)
  end

  new_coords
end
