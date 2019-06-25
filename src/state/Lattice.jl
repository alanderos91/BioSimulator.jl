struct Lattice{D,T,M,U}
  site::Vector{Site{D,T}}  # sorted by ID
  coord_order::Vector{Site{D,T}} # sort by coord
  neighbors::Vector{Vector{Int}}
  types::Vector{U}
end


function Lattice(coord::Matrix, types::Vector;
    nbhood = VonNeumann(),
    type_list = unique(types)
  )
  dimension = size(coord, 1)
  unique_types = sort!(type_list)
  number_types = length(unique_types)
  number_neighbors = capacity(nbhood, dimension)

  labels = OrderedDict(unique_types[i] => i+1 for i in eachindex(unique_types))
  site = [Site(i, State(labels[types[i]]), tuple(coord[:,i]...)) for i in eachindex(types)]

  site_by_coord = sort(site, by = coordinates)

  neighbors = [sizehint!(Int[], number_neighbors) for i in eachindex(site)]

  mapping = [(label => l) for (label, l) in labels]
  T = eltype(coord)
  U = eltype(mapping)

  return Lattice{dimension,T,typeof(nbhood),U}(site, site_by_coord, neighbors, mapping)
end

##### accessors

dimension(::Lattice{D}) where D = D
topology(::Lattice{D,T,M}) where {D,T,M} = M()

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

# this is a dumb hack to make Lattice compatible with SamplePath
Base.size(lattice::Lattice) = size(lattice.types)

function __simple_copy(lattice)
  coord = coordinates.(lattice.site)
  tcode = get_ptype.(lattice.site)

  idx = findall(!isequal(1), tcode)

  return coord[idx], tcode[idx]
end

##### query API

function istracked(lattice::Lattice, coord)
  idx = searchsorted(lattice.coord_order, coord, by = coordinates, lt = <)

  return !isempty(idx)
end

# get the site living at the given coordinates
function get_site(lattice::Lattice, coord)
  idx = searchsorted(lattice.coord_order, coord, by = coordinates, lt = <)

  lattice.coord_order[idx[1]]
end

# get the site with the given id
function get_site(lattice::Lattice, id::Int)
  # idx = searchsorted(lattice.site, id, by = label)
  #
  # # if idx is empty, site was not found
  #
  # # if idx has length > 1, we have a duplicate
  #
  # # if idx has length = 1, return the object
  # lattice.site[idx[1]]
  lattice.site[id]
end

function neighborhood(lattice::Lattice, x::Site)
  lattice.neighbors[label(x)]
end

##### update API

function spawn_new_site(lattice::Lattice, coord)
  id = number_sites(lattice) + 1

  return Site(id, State(), coord)
end

# this function assumes new_site is not tracked!
function add_site!(lattice::Lattice, new_site::Site)
  # add neighborhood for site
  nbmax = capacity(topology(lattice), dimension(lattice))
  push!(lattice.neighbors, sizehint!(Int[], nbmax))

  # new sites are always added to the end of the list sorted by IDs
  push!(lattice.site, new_site)

  # add site to list sorted by coord
  idx = searchsortedfirst(lattice.coord_order, new_site, by = coordinates)
  insert!(lattice.coord_order, idx[1], new_site)

  new_site
end

function add_neighbor!(lattice::Lattice{D}, x::Site{D}, y::Site{D}) where D
  nb = neighborhood(lattice, x)

  nbtype = topology(lattice)
  neighbor_count  = length(nb)
  nbhood_capacity = capacity(nbtype, D)

  if label(y) in nb
    error("$(y) is already a neighbor of $(x)")
  end

  if neighbor_count ≥ nbhood_capacity
    msg = """
    Neighborhood of $(x) at capacity ($(nbhood_capacity)).
    Failed to add $(y) as a neighbor.
    """
    nb_y = neighborhood(lattice, y)

    @debug """
    Neighbor data:

    x: $(label(x)) / $(get_ptype(x))|$(get_neighbor_class(x)) @ $(x)

    y: $(label(y)) / $(get_ptype(y))|$(get_neighbor_class(y)) @ $(y)

    Neighborhood:

    x:  $(nb)

    y:  $(nb_y)

    Sites:

    x:  $(lattice.site[nb])
    y:  $(lattice.site[nb_y])
    """
    throw(ErrorException(msg))
  end

  # store neighbors in sorted fashion
  i = label(y)
  idx = searchsortedfirst(nb, i)
  insert!(nb, idx, i)
end

function rmv_neighbor!(lattice::Lattice, x::Site, y::Site)
  i = label(x)
  j = label(y)

  nb = neighborhood(lattice, x)
  idx = findfirst(isequal(j), nb)

  if idx isa Nothing
    msg = """
    $(y) is not adjacent to $(x) under $(topology(lattice)) structure.
    Failed to remove neighbor.
    """
    throw(ErrorException(msg))
  end

  # delete item, which preserves ordering
  deleteat!(nb, idx)
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
    myseries = filter(x -> get_ptype(x) == last(t), site)

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
    s = filter(x -> get_ptype(x) == last(t), site)

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

##### building neighborhoods

function build_local_neighborhoods!(nbtype::VonNeumann, lattice::Lattice{1})
  sites = lattice.coord_order

  sort_by_x_1D!(sites)
  sweep_neighbors!(nbtype, lattice, sites)

  return nothing
end

function build_neighborhoods!(nbtype::VonNeumann, lattice::Lattice{2})
  sites = lattice.coord_order

  sort_by_x_2D!(sites)
  sweep_neighbors!(nbtype, lattice, sites)

  sort_by_y_2D!(sites)
  sweep_neighbors!(nbtype, lattice, sites)

  return nothing
end

function build_neighborhoods!(nbtype::VonNeumann, lattice::Lattice{3})
  sites = lattice.coord_order

  sort_by_x_3D!(sites)
  sweep_neighbors!(nbtype, lattice, sites)

  sort_by_y_3D!(sites)
  sweep_neighbors!(nbtype, lattice, sites)

  sort_by_z_3D!(sites)
  sweep_neighbors!(nbtype, lattice, sites)

  return nothing
end

function sweep_neighbors!(nbtype::VonNeumann, lattice, sites)
  for i = 2:length(sites)
    x = sites[i-1]
    y = sites[i]
    d = distance(nbtype, x, y)
    if d == 1 && !(label(y) in neighborhood(lattice, x))
      add_neighbor!(lattice, x, y)
      add_neighbor!(lattice, y, x)
    end
  end
  return sites
end

function build_neighborhoods!(nbtype::Hexagonal, lattice::Lattice{2})
  sites = lattice.coord_order

  sweep_neighbors!(nbtype, lattice, sites)
end

function sweep_neighbors!(nbtype::Hexagonal, lattice, sites)
  n = length(sites)

  for i in 1:n
    for j in 1:i-1 # j < i
      x = sites[i]
      y = sites[j]
      d = distance(nbtype, x, y)

      if d == 1 && !(label(y) in neighborhood(lattice, x))
        add_neighbor!(lattice, x, y)
        add_neighbor!(lattice, y, x)
      end
    end
  end
  return sites
end

### TODO: 3D
