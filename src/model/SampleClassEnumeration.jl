struct SampleClassEnumeration{D,M}
  class::Vector{Vector{Int}}
  composition::Vector{Vector{Int}}
  reactant_to_class::OrderedDict{Tuple{Int,Int},Int}
  pair_to_classes::OrderedDict{Tuple{Int,Int},Vector{Int}}
  isactive::Vector{Bool}
  dummy_composition::Vector{Int}
end

function SampleClassEnumeration{D,M}(composition, reactant_to_class, pair_to_classes, isactive) where {D,M}
  number_classes = reactant_to_class |> values |> maximum

  class = [Int[] for s in 1:number_classes]
  dummy_composition = zero(composition[1])

  SampleClassEnumeration{D,M}(
    class,
    composition,
    reactant_to_class,
    pair_to_classes,
    isactive,
    dummy_composition
  )
end

function build_reactant_pairs(reactions)
  # so, we need to build up the number of sample classes while forming the reaction channels.
  # this involves counting the number of unique reactants, singly for on-site, in pairs for pairwise.
  # ordering is important for the pairs too (though we can fix this later if we want...)
  pairs = Tuple{Int,Int}[]

  for reaction in reactions
    # gather up the pairs of reactants:
    push!(pairs, (reaction.input1, reaction.input2) )
  end

  # this returns the unique pairs of reactant species.  On-site reactions are denoted by the second reactant being `0`.
  pairs = unique(pairs)
end

function map_reactant_to_class(pairs, nbmax)
  # Now that we've got the pairs, we need to form the sample indices.
  # This is equivalent to enumerating the different
  # unique pairs of reactants, splitting the pairwise reactants into 2*d groups (1 for each possible number of
  # adjacent product particles) and leaving one group for each on-site reactant.  Each reaction (old reaction)
  # points to exactly one of these groups, the group that it samples from.

  # we need a way to turn a reactant pair (plus number of adjacent members of the second type) into a sample class
  # so we turn to Dictionaries.

  # probably ought to sizehint this too...
  # and add type restrictions...
  reactant_to_class = OrderedDict{Tuple{Int,Int},Int}()

  ### I'm going to do this the janky way, first thing that came to mind.  I'm sure there's a better way of doing this.

  counts = 1
  for pair in pairs
    # add the pair and the referenced sample class to the dictionary
    reactant_to_class[pair] = counts

    # check to see if the pair is pairwise or on-site
    # the second entry is `0` if the pair is on-site
    if pair[2] == 0
      # we only need to increase the count by 1, since only one
      # sampling index is needed
      counts += 1

      # if the pair is pairwise, need to increase the count by
      # nbmax, one for each possible number of neighbor reactants
    else
      counts += nbmax
    end
  end

  return reactant_to_class
end

function build_active_set(pairs, number_types)
  isactive = zeros(Bool, number_types + 1)

  for pair in pairs
    reactant = pair[1]
    isactive[reactant] = true
  end

  return isactive
end

function map_pair_to_classes(composition, reactant_to_class, isactive, number_types)
  pair_to_classes = OrderedDict{Tuple{Int,Int},Vector{Int}}()
  number_compositions = length(composition)

  for l in 1:(number_types + 1)
    # check if the species is active before doing anything
    if isactive[l]
      for k in 1:number_compositions
        ### get the composition associated with k, .
        comp = composition[k]

        # get an empty vector to collect the sample classes a particle of type l w/ neighborhood k takes part in
        sample_classes = Int[]

        # iterate over every reactant pair (l, j)
        iterate_reactant_pairs!(sample_classes, l, comp, reactant_to_class, number_types)

        # finally, put the sample_classes in the dictionary with key (l, k)
        pair_to_classes[(l, k)] = sample_classes
      end
    end
  end
  return pair_to_classes
end

function iterate_reactant_pairs!(sample_classes, l, comp, reactant_to_class, number_types)
  for j in 1:(number_types + 1)
    # get the number of neighbors of type j
    number_neighbors = comp[j]

    # if composition has type j...
    if number_neighbors != 0
      # get the sample class
      # return 0 if (l, j) is not a reacting pair
      s = get(reactant_to_class, (l, j), 0)

      # if (l,j) takes part in a reaction...
      if s != 0
        # push the appropriate sample class to the vector
        push!(sample_classes, s + number_neighbors - 1)
      end
    end
  end

  # finally, check to see if (l, 0) is a reactant pair
  # meaning l takes part in an on-site reaction
  s = get(reactant_to_class, (l, 0), 0)

  if s != 0
    # no need to correct s here
    push!(sample_classes, s)
  end

  return sample_classes
end

##### accessors

number_compositions(enum::SampleClassEnumeration) = length(enum.composition)

##### query API

# get the index assigned to a particular neighborhood class
function get_nbclass_index(enum::SampleClassEnumeration, comp)
  k = searchsortedfirst(enum.composition, comp, rev = true)

  if k > length(enum.composition)
    msg = """
    Composition $(comp) not found.
    """
    error(msg)
  end

  return k
end

# derive the neighbor class of a site with
# (i) nbmax - 1 open sites
# (ii) 1 neighbor of type i
function derive_nbclass(enum::SampleClassEnumeration{D,M}, i) where {D,M}
  # pull up the worker array
  comp = enum.dummy_composition

  # build the desired composition
  comp[1] = capacity(M(), D) - 1
  comp[i] = 1

  # search for the index; assumes composition is sorted somehow
  k = get_nbclass_index(enum, comp)

  # cleanup the mess
  fill!(comp, 0)

  return k
end

# get neighbor class after a i -> j type change,
# given that the site was in class k
function get_new_nbclass(enum::SampleClassEnumeration, k, i, j)
  # retrieve the list of compositions
  composition = enum.composition

  # retrieve the worker array
  comp = enum.dummy_composition

  # copy the composition we need to operate on
  copyto!(comp, composition[k])

  # generate the new composition
  comp[i] -= 1
  comp[j] += 1

  k_new = get_nbclass_index(enum, comp)

  # clean up
  fill!(comp, 0)

  return k_new
end

function add_class_member!(class, x::Site)
  i = label(x)
  idx = searchsortedfirst(class, i)
  insert!(class, idx, i)

  return nothing
end

function rmv_class_member!(class, x::Site)
  i = searchsortedfirst(class, label(x))
  deleteat!(class, i)

  return nothing
end

function sample_from_class(lattice, enum, s)
  return get_site(lattice, rand(enum.class[s]))
end

##### TODO Here be dragons.

# Update neighbor and sample classes for x, given that
# 	x changes from type i1 to j1
# 	y changes from type i2 to j2
function update_classes_particle!(enumeration, x, i1, j1, i2, j2)
  # unpack information
  k_old = get_neighbor_class(x)

  # 1. update neighborhood class if necessary
  if i2 != j2
    k_new = get_new_nbclass(enumeration, k_old, i2, j2)
    change_neighbor_class!(x, k_new)
  end

  # 2. update sample class
  # Case 1: i1 == j1 and i2 != j2; neighborhood structure changed, but not type
  if i1 == j1 && i2 != j2
    update_classes_neighbor_change!(enumeration, x, i1, k_old, i2, j2)
  end

  if i1 != j1
    # Case 2: i1 != j1 and i2 == j2; neighborhood structure for x did not change
    if i2 == j2
      # removals
      sample_class_removals!(enumeration, x, i1, k_old)

      # additions
      sample_class_additions!(enumeration, x, j1, k_old)
    else
    # Case 3: i1 != j1 and i2 != j2; neighborhood structure for x did change

      # removals
      sample_class_removals!(enumeration, x, i1, k_old)

      # additions
      sample_class_additions!(enumeration, x, j1, k_new)
    end
  end
  return nothing
end

function sample_class_removals!(enumeration, x, l, k)
  # unpack information
  isactive = enumeration.isactive
  pair_to_classes = enumeration.pair_to_classes
  class = enumeration.class

  if isactive[l]
    affected_classes = pair_to_classes[(l, k)]

    for s in affected_classes
      rmv_class_member!(class[s], x)
    end
  end

  return enumeration
end

function sample_class_additions!(enumeration, x, l, k)
  # unpack information
  isactive = enumeration.isactive
  pair_to_classes = enumeration.pair_to_classes
  class = enumeration.class

  if isactive[l]
    affected_classes = pair_to_classes[(l, k)]

    for s in affected_classes
      add_class_member!(class[s], x)
    end
  end

  return enumeration
end

function update_classes_neighbor_change!(enumeration, x, i1, k_old, i2, j2)
  isactive = enumeration.isactive
  reactant_to_class = enumeration.reactant_to_class
  class = enumeration.class
  composition = enumeration.composition

  if isactive[i1]
    # what class corresponds to the reactant pair (i1, i2)?
    s = get(reactant_to_class, (i1, i2), 0)

    # updates for losing particle of type i2
    if s != 0
      # get number of neighbors of type i2, the interacting particle
      number_neighbors = composition[k_old][i2]

      # index shifting
      s += number_neighbors - 1

      # remove from the sample class
      rmv_class_member!(class[s], x)

      if number_neighbors != 1
        # lost a particle of type i2, so add it to the previous class
        add_class_member!(class[s-1], x)
      end
    end

    # updates for gaining a particle of type j2
    s = get(reactant_to_class, (i1, j2), 0)

    if s != 0
      number_neighbors = composition[k_old][j2]

      s += number_neighbors

      add_class_member!(class[s], x)

      if number_neighbors != 0
        rmv_class_member!(class[s-1], x)
      end
    end
  end

  return nothing
end

# update neighbors due to x interacting with y, where x changes from i to j
function update_classes_neighbors!(lattice, enumeration, x, y, i, j)
  # N = lattice.N
  nbtype = topology(lattice)

  # check the coordinates around x under the given neighborhood type
  for coord in eachdir(nbtype, x)
    # if there is a particle at the coordinate...
    if istracked(lattice, coord)
      neighbor_site = get_site(lattice, coord)

      # if the neighbor site is not the one that interacted with x...
      if neighbor_site != y
        # change neighborhood classes
        k_old = get_neighbor_class(neighbor_site)
        k_new = get_new_nbclass(enumeration, k_old, i, j)

        change_neighbor_class!(neighbor_site, k_new)

        # change sample classes
        l = get_ptype(neighbor_site)
        update_classes_neighbor_change!(enumeration, neighbor_site, l, k_old, i, j)

        # 3. update counts only if neighbor_site did not reaction with x
        # N[k_old, l] -= 1
        # N[k_new, l] += 1
      end
    # otherwise, there is an open site at the coordinate
    else
      neighbor_site = spawn_new_site(lattice, coord)
      l = get_ptype(neighbor_site)

      # set neighborhood class and assign sample classes
      k = derive_nbclass(enumeration, j)
      change_neighbor_class!(neighbor_site, k)
      sample_class_additions!(enumeration, neighbor_site, l, k)

      # add the site
      add_site!(lattice, neighbor_site)

      # update the neighborhood
      for coord2 in eachdir(nbtype, neighbor_site)
        if istracked(lattice, coord2)
          z = get_site(lattice, coord2)

          is_neighbor = label(z) in neighborhood(lattice, neighbor_site)

          if !is_neighbor
            add_neighbor!(lattice, neighbor_site, z)
            add_neighbor!(lattice, z, neighbor_site)
          end
        end
      end
    end
  end

  return nothing
end

# we could make this use a binary search if neighbors
# are stored in a sorted fashion
function sample_neighbor(x, lattice, enum, l)
  # build an iterable object that returns the type of site for each neighbor
  site = lattice.site
  nb   = neighborhood(lattice, x)
  iter = (get_ptype(site[id]) for id in nb)

  # extract the number of type l neighbors
  k = get_neighbor_class(x)
  number_l = enum.composition[k][l]

  # generate a random index into nb
  ix = sample_neighbor(iter, number_l, l)

  if ix == 0
    msg = """
    $(get_ptype(x))|$(get_neighbor_class(x)) @ $(label(x)) / $(x) does not have a neighbor of type $(l).
      Adjacent types: $(map(get_ptype, lattice.site[nb]))
    """
    error(msg)
  end

  # return the chosen Site object
  get_site(lattice, nb[ix])
end

function sample_neighbor(iter, number_l, l)
  c = number_l * rand()
  s = 0 # tracks the number of type l neighbors we've checked
  j = 0 # tracks the random index we want to generate

  for val in iter
    s > c && break
    isequal(val, l) && (s += 1)
    j += 1
  end

  return j
end

# # TESTING
#
# number_trials = 10_000
# x = rand(0:9, 20)
# l = rand(0:9)
# number_l = count(isequal(l), x)
# iter = (x[j] for j in eachindex(x))
#
# # 1. we always return an index mapping to l
#
# for k in 1:number_trials
#   i = sample_neighbor(iter, number_l, l)
#   @assert x[i] == l
# end
#
# # 2. every index is mapping to l is touched
#
# idxs = Int[]
#
# for k in 1:number_trials
#    push!(idxs, sample_neighbor(iter, number_l, l))
# end
# @assert length(unique(idxs)) == number_l
#
# # 3. the distribution should be uniform
# histogram(idxs)
#
