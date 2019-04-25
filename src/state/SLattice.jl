# ptype == 1 is reserved for open sites!
struct SLattice{D,T,M} <: AbstractLattice{D,T,M}
  classes::Vector{SampleClass{Site{D,T,M}}}
  k2composition::OrderedDict{Int,Vector{Int}}
  composition2k::OrderedDict{Vector{Int},Int}
  label2site::OrderedDict{Int,Site{D,T,M}}
  coord2site::OrderedDict{NTuple{D,T},Site{D,T,M}}

  # new fields we need
  isactive::Vector{Bool}
  reac2sidx::OrderedDict{Tuple{Int,Int},Int}
  lk2sidx::OrderedDict{Tuple{Int,Int},Vector{Int}}

  # matrix of particle counts
  N :: Matrix{Int}

  # these fields ought to become type parameters in the future
  L :: Int
  K :: Int

  dummy_composition :: Vector{Int}

  # xlim :: Tuple{Int,Int}
  # ylim :: Tuple{Int,Int}
  # zlim :: Tuple{Int,Int}
  
  function SLattice{D,T,M}(sites::Vector{Site{D,T,M}}, compositions, L, K, reactions) where {D,T,M}
    # map an index to a composition
    k2composition = OrderedDict{Int,Vector{Int}}(i => compositions[i] for i in eachindex(compositions))

    # map a composition to an index
    composition2k = OrderedDict{Vector{Int},Int}(compositions[i] => i for i in eachindex(compositions))

    HINT = 1_000

    # map a label to a `Site`
    label2site = OrderedDict{Int,Site{D,T,M}}(label(site) => site for site in sites)
    sizehint!(label2site, HINT)

    # map a coordinate to a `Site`
    coord2site = OrderedDict{NTuple{D,T},Site{D,T,M}}(tuple(coordinates(site)...) => site for site in sites)
    sizehint!(coord2site, HINT)

    # add missing open sites
    m = length(sites)

    # build the reactant pairs
    pairs = build_reactant_pairs(reactions)

    # build the reactant pair to sample index dictionary
    nbmax = capacity(neighborhood(sites[1]))
    reac2sidx = build_pair_dictionary(pairs, nbmax)

    # build active set
    isactive = build_active_set(pairs, L)

    # build the (tuple, neighborhood class) to sample index dicitionary
    lk2sidx = build_lk2sidx_dict(compositions, reac2sidx, isactive, L)

    for i in 1:m
      site = sites[i]
      for coord in eachdir(site)
        if !haskey(coord2site, coord)
          id = length(coord2site) + 1
          neighbor_site = Site(id, State(), coord, M)

          coord2site[coord] = neighbor_site
          label2site[id] = neighbor_site
          push!(sites, neighbor_site)
        end
      end
    end

    # build local neighborhoods
    build_local_neighborhoods!(sites)

    # determine number of sample classes
    number_pairwise_reactants = count(reac_pair -> reac_pair[2] != 0, pairs)
    number_onsite_reactants   = count(reac_pair -> reac_pair[2] == 0, pairs)
    number_sample_classes     = 2 * D * number_pairwise_reactants + number_onsite_reactants
    
    # build the sample classes
    classes = [SampleClass{Site{D,T,M}}(s) for s in 1:number_sample_classes]

    composition = zeros(Int, L + 1)
    dummy_empty_list = Int[]
    for site in sites
      build_composition!(composition, site)

      l = get_ptype(site)
      k = composition2k[composition]

      # this will exclude empty sites from being added into a sample class
      sample_classes = get(lk2sidx, (l,k), dummy_empty_list)

      for s in sample_classes
        add_member!(classes[s], site)
      end

      # don't forget to set the correct neighborhood class
      change_neighbor_class!(site, k)
      
      # add_member!(classes[k], site)
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

    return new{D,T,M}(classes, k2composition, composition2k, label2site, coord2site, isactive, reac2sidx, lk2sidx, N, L, K, dummy_composition)
  end
end

# function barrier to deal with dimensionality
function SLattice(p::Vector{Tuple{NTuple{D,Int64},Int64}}, L, reactions) where D
  compositions = collect(Vector{Int64}, multiexponents(L + 1, 2 * D))
  sites = [Site(i, State(p[i][2]), p[i][1]) for i in eachindex(p)]
  K = length(compositions)

  return SLattice(sites, compositions, L, K, reactions)
end

# generalize to different neighborhood types
function SLattice(p::Vector{Tuple{NTuple{D,Int64},Int64}}, L, reactions, nbhood) where D
  compositions = collect(Vector{Int64}, multiexponents(L + 1, 2 * D))
  
  if nbhood == :vonneumann
    if D == 1
      M = VonNeumann1D{Int}
    elseif D == 2
      M = VonNeumann2D{Int}
    else
      M = VonNeumann3D{Int}
    end
  elseif nbhood == :hexagonal
    if D == 2
      M = Hexagonal2D{Int}
    end
  end

  sites = [Site(i, State(p[i][2]), p[i][1], M) for i in eachindex(p)]
  K = length(compositions)

  return SLattice(sites, compositions, L, K, reactions)
end

# where does this fall?
SLattice(sites::Vector{Site{D,T,M}}, compositions, L, K, reactions) where {D,T,M} = SLattice{D,T,M}(sites, compositions, L, K, reactions)

# this is how users build a lattice
function SLattice(point::Matrix{Int}, ptype::Vector{Int}, reactions; number_types = nothing, nbhood = :vonneumann)
  # if the user does not specify the number of types, then we assume all the
  # types occur inside the ptype list
  if number_types == nothing
    L = length(unique(ptype))
  else
    L = number_types :: Int
  end
  
  p = [(tuple(point[:, i]...), ptype[i] + 1) for i in eachindex(ptype)]

  return SLattice(p, L, reactions, nbhood)
end

# update api #

# sample a site from sample class s
sample_from_class(lattice::SLattice, s) = sample(lattice.classes[s])

function add_site!(lattice::SLattice, x)
  k = get_neighbor_class(x)
  
  # add to neighborhood class; don't add open sites for now
  if get_ptype(x) > 1
    add_member!(lattice.classes[k], x)
  end

  # add to label2site map
  lattice.label2site[label(x)] = x

  # add to coord2site map
  lattice.coord2site[tuple(coordinates(x)...)] = x

  return nothing
end


# Update neighbor and sample classes for x, given that
# 	x changes from type i1 to j1
# 	y changes from type i2 to j2
function update_classes_particle!(lattice::SLattice, x, i1, j1, i2, j2)
  # unpack information
  k_old = get_neighbor_class(x)

  # 1. update neighborhood class if necessary
  if i2 != j2
    k_new = get_new_class(lattice, k_old, i2, j2)
    change_neighbor_class!(x, k_new)
  end

  # 2. update sample class
  # Case 1: i1 == j1 and i2 != j2; neighborhood structure changed, but not type
  # use reac2sidx to get sample classes
  if i1 == j1 && i2 != j2
    update_classes_neighbor_change!(lattice, x, i1, k_old, i2, j2)
  end

  if i1 != j1
    # Case 2: i1 != j1 and i2 == j2; neighborhood structure for x did not change
    if i2 == j2
      # removals
      sample_class_removals!(lattice, x, i1, k_old)

      # additions
      sample_class_additions!(lattice, x, j1, k_old)
    else
    # Case 3: i1 != j1 and i2 != j2; neighborhood structure for x did change

      # removals
      sample_class_removals!(lattice, x, i1, k_old)

      # additions
      sample_class_additions!(lattice, x, j1, k_new)
    end
  end
  return nothing
end

function sample_class_removals!(lattice::SLattice, x, l, k)
  if lattice.isactive[l]
    classes = lattice.lk2sidx[(l, k)]

    for s in classes
      rmv_member!(lattice.classes[s], x)
    end
  end
  return lattice
end

function sample_class_additions!(lattice::SLattice, x, l, k)
  if lattice.isactive[l]
    classes = lattice.lk2sidx[(l, k)]

    for s in classes
      add_member!(lattice.classes[s], x)
    end
  end
  return lattice
end

function update_classes_neighbor_change!(lattice::SLattice, x, i1, k_old, i2, j2)
  if lattice.isactive[i1]
    s = get(lattice.reac2sidx, (i1, i2), 0)

    # updates for losing particle of type i2
    if s != 0
      # get number of neighbors of type i2, the interacting particle
      number_neighbors = lattice.k2composition[k_old][i2]

      # index shifting
      s += number_neighbors - 1

      # remove from the sample class
      rmv_member!(lattice.classes[s], x)

      if number_neighbors != 1
        # lost a particle of type i2, so add it to the previous class
        add_member!(lattice.classes[s - 1], x)
      end
    end

    # updates for gaining a particle of type j2
    s = get(lattice.reac2sidx, (i1, j2), 0)

    if s != 0
      number_neighbors = lattice.k2composition[k_old][j2]

      s += number_neighbors

      add_member!(lattice.classes[s], x)

      if number_neighbors != 0
        rmv_member!(lattice.classes[s - 1], x)
      end
    end
  end
  return nothing
end

# update neighbors due to x interacting with y, where x changes from i to j
function update_classes_neighbors!(lattice::SLattice, x, y, i, j)
  N = lattice.N

  for coord in eachdir(x)
    if istracked(lattice, coord)
      neighbor_site = get_site(lattice, coord)
      if neighbor_site != y
        # 1. change neighborhood classes
        k_old = get_neighbor_class(neighbor_site)
        k_new = get_new_class(lattice, k_old, i, j)

        # println()
        # println("tried to update $(neighbor_site)")
        # println("$(k_old) ----> $(k_new)")
        # println()

        change_neighbor_class!(neighbor_site, k_new)

        # 2. change sample classes
        l = get_ptype(neighbor_site)
        update_classes_neighbor_change!(lattice, neighbor_site, l, k_old, i, j)

        # 3. update counts only if neighbor_site did not reaction with x
        N[k_old, l] -= 1
        N[k_new, l] += 1
      end
    else
      neighbor_site = spawn_new_site(lattice, coord)
      l = get_ptype(neighbor_site)

      # set particle class
      k = derive_class(lattice, j)
      change_neighbor_class!(neighbor_site, k)
      sample_class_additions!(lattice, neighbor_site, l, k)
      
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
      N[k, l] += 1

      # add the site
      add_site!(lattice, neighbor_site)
    end
  end

  return nothing
end


