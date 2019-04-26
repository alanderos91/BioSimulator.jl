# type declaration #
abstract type AbstractLattice{D,T,M} end

# check boundary #

@inline function inboundary(coord :: NTuple{1,T}, xlim :: NTuple{1,T}, ylim :: NTuple{1,T}) where T
  return xlim[1] <= coord[1] <= xlim[2]
end

@inline function inboundary(coord :: NTuple{2,T}, xlim :: NTuple{2,T}, ylim :: NTuple{2,T}) where T
  return (xlim[1] <= coord[1] <= xlim[2]) && (ylim[1] <= coord[2] <= ylim[2])
end

@inline function inboundary(coord :: NTuple{3,T}, xlim :: NTuple{3,T}, ylim :: NTuple{3,T}) where T
  return (xlim[1] <= coord[1] <= xlim[2]) && (ylim[1] <= coord[2] <= ylim[2]) && (zlim[1] <= coord[3] <= zlim[2])
end

# accessors #
number_of_types(lattice::AbstractLattice) = lattice.L
number_of_compositions(lattice::AbstractLattice) = lattice.K

# update API #

function spawn_new_site(lattice::AbstractLattice{D,T,M}, coord) where {D,T,M}
	label = length(lattice.label2site) + 1
	site = Site(label, State(), coord, M)

	return site
end

# get the site living at the given coordinate
get_site(lattice::AbstractLattice, coordinate) = lattice.coord2site[coordinate]

# is the site in our dictionary?
istracked(lattice::AbstractLattice, coordinate) = haskey(lattice.coord2site, coordinate)

# this computes the class of a composition with C_0 = 2 * D - 1, C_i = 1
function derive_class(lattice::AbstractLattice{D,T,M}, i) where {D,T,M}
  # we can improve on this using StaticArrays...
  composition = lattice.dummy_composition

  composition[1] = capacity(M) - 1
  composition[i] = 1

  k = lattice.composition2k[composition]
  
  # clean up
  fill!(composition, 0)

  return k
end

# was in class k, lost neighbor of type i that became type j
function get_new_class(lattice::AbstractLattice, k, i, j)
  # retrieve the dummy composition
  composition = lattice.dummy_composition

  # copy the composition we need to operate on
  copyto!(composition, lattice.k2composition[k])

  # println("old composition = ", composition)
  # generate the new composition
  composition[i] -= 1;
  composition[j] += 1
  # println("new composition = ", composition)

  m = lattice.composition2k[composition]

  # clean up
  fill!(composition, 0)

  return m
end

# dirty hacks

# length(lattice::AbstractLattice) = length(lattice.N)
Base.copy(lattice::AbstractLattice) = deepcopy(lattice)