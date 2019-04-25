# ptype == 1 is reserved for open sites!
struct NLattice{D,T,M} <: AbstractLattice{D,T,M}
  classes::Vector{NeighborhoodClass}
  k2composition::Dict{Int,Vector{Int}}
  composition2k::Dict{Vector{Int},Int}
  label2site::Dict{Int,Site{D,T,M}}
  coord2site::Dict{NTuple{D,T},Site{D,T,M}}
  
  # matrix of particle counts
  N :: Matrix{Int}
  
  # these fields ought to become type parameters in the future
  L :: Int
  K :: Int
  
  dummy_composition :: Vector{Int}
  
  # xlim :: Tuple{Int,Int}
  # ylim :: Tuple{Int,Int}
  # zlim :: Tuple{Int,Int}
  
  function NLattice{D,T,M}(sites::Vector{Site{D,T,M}}, compositions, L, K) where {D,T,M}
    # map an index to a composition
    k2composition = Dict{Int,Vector{Int}}(i => compositions[i] for i in eachindex(compositions))
    
    # map a composition to an index
    composition2k = Dict{Vector{Int},Int}(compositions[i] => i for i in eachindex(compositions))
    
    HINT = 1_000
    
    # map a label to a `Site`
    label2site = Dict{Int,Site{D,T,M}}(label(site) => site for site in sites)
    sizehint!(label2site, HINT)
    
    # map a coordinate to a `Site`
    coord2site = Dict{NTuple{D,T},Site{D,T,M}}(tuple(coordinates(site)...) => site for site in sites)
    sizehint!(coord2site, HINT)
    
    # add missing open sites
    m = length(sites)
    
    for i in 1:m
      site = sites[i]
      for coord in eachdir(site)
        if !haskey(coord2site, coord)
          id = length(coord2site) + 1
          neighbor_site = Site(id, State(), coord)
          
          coord2site[coord] = neighbor_site
          label2site[id] = neighbor_site
          push!(sites, neighbor_site)
        end
      end
    end
    
    # build local neighborhoods
    build_local_neighborhoods!(sites)
    
    sort!(sites, by = label)
    
    # build the neighborhood classes
    classes = [NeighborhoodClass{Site{D,T,M}}(k) for k in eachindex(compositions)]
    
    composition = zeros(Int, L + 1)
    
    for site in sites
      build_composition!(composition, site)
      k = composition2k[composition]
      
      change_neighbor_class!(site, k)
      
      add_member!(classes[k], site)
    end
    
    # build the count vector
    N = zeros(Int, K, L + 1)
    
    for site in sites
      k = get_neighbor_class(site)
      l = get_ptype(site)
      
      N[k, l] += 1
    end
    
    # xlim = ?
    # ylim = ?
    # zlim = ?
    
    # return new(classes, k2composition, composition2k, label2site, coord2site, N, L, K, xlim, ylim, zlim)
    
    # make a dummy vector for doing manipulations with compositions
    dummy_composition = zeros(Int, L + 1)
    
    return new{D,T,M}(classes, k2composition, composition2k, label2site, coord2site, N, L, K, dummy_composition)
  end
end

# function barrier to deal with dimensionality
function NLattice(p::Vector{Tuple{NTuple{D,Int64},Int64}}, L) where D
  compositions = collect(Vector{Int64}, multiexponents(L + 1, 2 * D))
  sites = [Site(i, State(p[i][2]), p[i][1]) for i in eachindex(p)]
  K = length(compositions)
  
  return NLattice(sites, compositions, L, K)
end

NLattice(sites::Vector{Site{D,T,M}}, compositions, L, K) where {D,T,M} = NLattice{D,T,M}(sites, compositions, L, K)

# this is how users build a lattice
function NLattice(point::Matrix{Int}, ptype::Vector{Int}; number_types = nothing)
  # if the user does not specify the number of types, then we assume all the
  # types occur inside the ptype list
  if number_types == nothing
    L = length(unique(ptype))
  else
    L = number_types :: Int
  end
  
  p = [(tuple(point[:, i]...), ptype[i] + 1) for i in eachindex(ptype)]
  
  return NLattice(p, L)
end

# TODO

# update api #

# sample a site of type l in neighborhood class k
function sample_center(lattice::NLattice, l, k)
  sample_center(lattice.classes[k], l, lattice.N[k, l])
end

function change_class!(lattice::NLattice, x::Site, k_old, k_new)
  # grab the two neighborhood classes
  nc_old = lattice.classes[k_old]
  nc_new = lattice.classes[k_new]
  
  # remove x from its old neighborhood class, then add it to the new one
  rmv_member!(nc_old, x)
  add_member!(nc_new, x)
  
  return nothing
end

function add_site!(lattice::NLattice, x)
  k = get_neighbor_class(x)
  
  # add to neighborhood class
  add_member!(lattice.classes[k], x)
  
  # add to label2site map
  lattice.label2site[label(x)] = x
  
  # add to coord2site map
  lattice.coord2site[tuple(coordinates(x)...)] = x
  
  return nothing
end

function update_neighbor_classes!(lattice::NLattice, x::Site, y::Site, i, j)
  N = lattice.N
  
  for coord in eachdir(x)
    if istracked(lattice, coord)
      neighbor_site = get_site(lattice, coord)
      
      k_old = get_neighbor_class(neighbor_site)
      k_new = get_new_class(lattice, k_old, i, j)
      
      # println()
      # println("tried to update $(neighbor_site)")
      # println("$(k_old) ----> $(k_new)")
      # println()
      
      change_neighbor_class!(neighbor_site, k_new)
      change_class!(lattice, neighbor_site, k_old, k_new)
      
      # update counts only if neighbor_site did not interact with x
      if neighbor_site != y
        N[k_old, get_ptype(neighbor_site)] -= 1
        N[k_new, get_ptype(neighbor_site)] += 1
      end
    else
      neighbor_site = spawn_new_site(lattice, coord)
      
      # set particle class
      k = derive_class(lattice, j)
      change_neighbor_class!(neighbor_site, k)
      
      # local neighborhood
      # add_neighbor!(x, neighbor_site)
      # add_neighbor!(neighbor_site, x)
      for coord2 in eachdir(neighbor_site)
        if istracked(lattice, coord2)
          z = get_site(lattice, coord2)
          if z ∉ neighborhood(neighbor_site)
            add_neighbor!(neighbor_site, z)
            add_neighbor!(z, neighbor_site)
          end
        end
      end 
      
      # update counts
      N[k, get_ptype(neighbor_site)] += 1
      
      # add the site
      add_site!(lattice, neighbor_site)
    end
  end
  
  return nothing
end
