struct Lattice{D,T,M,U}
  site::Vector{Site{D,T}}
  neighbors::Vector{Vector{Int}}
  types::Vector{U}
end

function Lattice(coord::Matrix{T}, types::Vector{U}; nbhood = VonNeumann()) where {T,U}
  dimension = size(coord, 1)
  number_types = length(unique(types))
  number_neighbors = capacity(nbhood, dimension)

  site = [Site(i, State(types[i]), tuple(coord[:,i]...)) for i in eachindex(types)]

  neighbors = [sizehint!(Int[], number_neighbors) for i in eachindex(site)]

  return Lattice{dimension,T,typeof(nbhood),U}(site, neighbors, types)
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
