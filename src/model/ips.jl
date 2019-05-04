struct IPSReactionStruct
  pairwise :: Bool
  input1   :: Int # center
  input2   :: Int # adjacent
  output1  :: Int # center
  output2  :: Int # adjacent
  class    :: Int
  rate     :: Float64
  # formula    :: Expr
  # parameter  :: Symbol
  # paramidx   :: Int

  function IPSReactionStruct(pairwise, i1, i2, o1, o2, class, rate)
    if has_invalid_type(pairwise, i1, i2, o1, o2)
      throw(ErrorException("""
      Particle types must be positive integers. Non-pairwise reactions must have the adjacent types be zero.

      pairwise? $pairwise
      center:   $i1 -> $o1
      adjacent: $i2 -> $o2
      """))
    end

    if has_same_types(pairwise, i1, i2, o1, o2)
      throw(ErrorException("""
      At least one particle type must be different.

      pairwise? $pairwise
      center:   $i1 -> $o1
      adjacent: $i2 -> $o2
      """))
    end

    if has_invalid_class(class)
      throw(ErrorException("Assigned class must be a positive integer; class = $class"))
    end

    if has_invalid_rate(rate)
      throw(ErrorException("Reaction rate must be a non-negative real number; rate = $rate"))
    end

    return new(pairwise, i1, i2, o1, o2, class, rate)
  end
end

@inline function has_invalid_type(pairwise, i1, i2, o1, o2)
  if pairwise
    # all types must be positive integers
    is_invalid = (i1 < 1) || (i2 < 1) || (o1 < 1) || (o2 < 1)
  else
    # center types must be positve
    # adjacent type must be zero
    is_invalid = (i1 < 1) || (o1 < 1) || (i2 != 0) || (o2 != 0)
  end

  return is_invalid
end

@inline function has_same_types(pairwise, i1, i2, o1, o2)
  if pairwise
    is_invalid = (i1 == i2) && (i2 == o1) && (o1 == o2)
  else
    is_invalid = (i1 == o1)
  end

  return is_invalid
end

@inline function has_invalid_class(class)
  # class must be a positive integer
  class < 1
end

@inline function has_invalid_rate(rate)
  # rate must be a non-negative real number
  rate < 0
end

ispairwise(x::IPSReactionStruct) = x.pairwise

@inline get_reactant_types(x::IPSReactionStruct) = (x.input1, x.input2)
@inline get_product_types(x::IPSReactionStruct) = (x.output1, x.output2)

##### pretty printing #####

function _print_formula(x::IPSReactionStruct)
  ex = formula(x)
  return string(ex.args[1], " ", "→", " ", ex.args[2])
end

struct InteractingParticleSystem{DG}
  reactions::Vector{IPSReactionStruct}
  # rxn_rates::Vector{T}
  dep_graph::DG
  # spc_graph::DG
  # rxn_graph::DG
end

function InteractingParticleSystem(reactions::Vector{IPSReactionStruct})
  num_reactions = length(reactions)

  dep_graph = rxnrxn_depgraph(DGLazy(), reactions)

  # return InteractingParticleSystem(reactions, rxn_rates, dep_graph, spc_graph, rxn_graph)
  return InteractingParticleSystem(reactions, dep_graph)
end

##### execute_jump! #####

execute_jump!(state, model::InteractingParticleSystem, j) = execute_jump!(state, model.reactions[j])

function execute_jump!(lattice::AbstractLattice, reaction::IPSReactionStruct)
  if ispairwise(reaction)
    # println("pairwise")
    execute_pairwise!(lattice, reaction)
  else
    # println("onsite")
    execute_onsite!(lattice, reaction)
  end

  return nothing
end

##### SLattice subroutines #####

function execute_pairwise!(lattice::SLattice, reaction)
  # unpack reaction information
  i1, i2 = get_reactant_types(reaction)
  j1, j2 = get_product_types(reaction)

  s = reaction.class

  # unpack counts matrix
  N = lattice.N

  # println("reaction = $(reaction)")

  # sample a particle x of type i1 from the given class k
  # println("attempt to sample $(i1) from class $(k)")
  x = sample_from_class(lattice, s)

  # sample a particle y of type i2 from x's neighborhood
  # println("attempt to sample $(i2) from neighborhood of $(x)")
  # println()
  y = sample(neighborhood(x), i2, lattice.k2composition[get_neighbor_class(x)][i2])

  #### this where we start ####

  # we use up the particles (k1, i1) and (k2, i2)

  N[get_neighbor_class(x), get_ptype(x)] -= 1
  N[get_neighbor_class(y), get_ptype(y)] -= 1

  # 1. change the types for x and y
  i1 != j1 && transform!(x, j1)
  i2 != j2 && transform!(y, j2)

  # 2. update neighborhood and sample classes for x and y
  update_classes_particle!(lattice, x, i1, j1, i2, j2)
  update_classes_particle!(lattice, y, i2, j2, i1, j1)

  # 3. update the remaining particles surrounding x and y
  i1 != j1 && update_classes_neighbors!(lattice, x, y, i1, j1)
  i2 != j2 && update_classes_neighbors!(lattice, y, x, i2, j2)

  # we produced new particles so update the counts
  N[get_neighbor_class(x), get_ptype(x)] += 1
  N[get_neighbor_class(y), get_ptype(y)] += 1

  return nothing
end

function execute_onsite!(lattice::SLattice, reaction)
  # unpack reaction information
  i, _ = get_reactant_types(reaction)
  j, _ = get_product_types(reaction)
  s = reaction.class

  # unpack counts matrix
  N = lattice.N

  # sample a particle x of type i from the given class k
  # println("attempt to sample $(i) from class $(k)")
  x = sample_from_class(lattice, s)

  # we use up the x particle
  N[get_neighbor_class(x), get_ptype(x)] -= 1

  # 1. change the type for x
  transform!(x, j)

  # 2. update neighborhood and sample classes for x and y
  update_classes_particle!(lattice, x, i, j, 0, 0)

  # 3. update the remaining particles surrounding x and y
  update_classes_neighbors!(lattice, x, x, i, j)

  # update the counts
  N[get_neighbor_class(x), get_ptype(x)] += 1

  return nothing
end

function rate(reaction::IPSReactionStruct, lattice::SLattice)
  # grab reactant indices
  l1, l2 = get_reactant_types(reaction)
  s = reaction.class

  # grab per particle rate
  per_particle_rate = reaction.rate

  # grab the count for the center particle
  N_s = length(lattice.classes[s].members) # need to fix this

  return N_s * per_particle_rate
end

##### NLattice subroutines #####

function execute_pairwise!(lattice::NLattice, reaction)
  # unpack reaction information
  i1, i2 = get_reactant_types(reaction)
  j1, j2 = get_product_types(reaction)
  k = reaction.class

  # unpack counts matrix
  N = lattice.N

  # println("reaction = $(reaction)")

  # sample a particle x of type i1 from the given class k
  # composition = lattice.k2composition[k]
  # println("  attempt to sample $(i1) from neighborhood class $(k) / $(composition)")
  x = sample_center(lattice, i1, k)
  # println("    got $(x)")
  # sample a particle y of type i2 from x's neighborhood
  # kx = get_neighbor_class(x)
  # composition = lattice.k2composition[kx]
  # println("  attempt to sample $(i2) from neighborhood class $(kx) / $(composition)")

  y = sample(neighborhood(x), i2, lattice.k2composition[k][i2])

  # println("    got $(y)")
  # println()

  # we use up the particles (k1, i1) and (k2, i2)
  N[get_neighbor_class(x), get_ptype(x)] -= 1
  N[get_neighbor_class(y), get_ptype(y)] -= 1

  # if the type of center particle changed...
  if i1 != j1
    # update the class of its neighbors
    # lost type i1 that became type j1
    update_neighbor_classes!(lattice, x, y, i1, j1)
  end

  # if the type of the neighboring particle changed...
  if i2 != j2
    # update the class of its neighbors
    # lost type i2 that became type j2
    update_neighbor_classes!(lattice, y, x, i2, j2)
  end

  # change states
  transform!(x, j1)
  transform!(y, j2)

  # we produced new particles so update the counts
  N[get_neighbor_class(x), get_ptype(x)] += 1
  N[get_neighbor_class(y), get_ptype(y)] += 1

  return nothing
end

function execute_onsite!(lattice::NLattice, reaction)
  # unpack reaction information
  i, _ = get_reactant_types(reaction)
  j, _ = get_product_types(reaction)
  k = reaction.class

  # unpack counts matrix
  N = lattice.N

  # sample a particle x of type i from the given class k
  # println("attempt to sample $(i) from class $(k)")
  x = sample_center(lattice, i, k)

  # we use up the x particle
  N[get_neighbor_class(x), get_ptype(x)] -= 1

  # we might need this later if we model some complicated internal state
  # # if the type of center particle changed...
  # if i != j
  #   # update the class of its neighbors
  #   # lost type i1 that became type j1
  #   update_neighbor_classes!(lattice, x, i, j)
  # end

  # for now, we assume i != j
  update_neighbor_classes!(lattice, x, x, i, j)

  # change states
  transform!(x, j)

  # update the counts
  N[get_neighbor_class(x), get_ptype(x)] += 1

  return nothing
end

function rate(reaction::IPSReactionStruct, lattice::NLattice)
  # grab reactant indices
  l1, l2 = get_reactant_types(reaction)
  k = reaction.class # neighbor class!

  # grab per particle rate
  per_particle_rate = reaction.rate

  # grab the count for the center particle
  N_kl = lattice.N[k, l1]

  return N_kl * per_particle_rate
end

##### convenience wrapper
rate(model::InteractingParticleSystem, lattice, j) = rate(model.reactions[j], lattice)

function build_reactant_pairs(reactions)
  # so, we need to build up the number of sample classes while forming the reaction channels.
  # this involves counting the number of unique reactants, singly for on-site, in pairs for pairwise.
  # ordering is important for the pairs too (though we can fix this later if we want...)
  pairs = Tuple{Int,Int}[]

  # sizehint `pairs`?
  for reaction in reactions
    # gather up the pairs of reactants:
    push!(pairs, (reaction.input1, reaction.input2) )
  end

  # this returns the unique pairs of reactant species.  On-site reactions are denoted by the second reactant being `0`.
  pairs = unique(pairs)
end

function build_pair_dictionary(pairs, nbmax)
  # now that we've got the pairs, we need to form the sample indices.  This is equivalent to enumerating the different
  # unique pairs of reactants, splitting the pairwise reactants into 2*d groups (1 for each possible number of
  # adjacent product particles) and leaving one group for each on-site reactant.  Each reaction (old reaction)
  # points to exactly one of these groups, the group that it samples from.


  # we need a way to turn a reactant pair (plus number of adjacent members of the second type) into a sample class
  # so we turn to Dictionaries.

  # probably ought to sizehint this too...
  # and add type restrictions...
  reac2sidx = OrderedDict{Tuple{Int,Int},Int}()


  ### I'm going to do this the janky way, first thing that came to mind.  I'm sure there's a better way of doing this.

  counts = 1
  for pair in pairs
    # add the pair and the referenced sample class to the dictionary
    reac2sidx[pair] = counts

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

  return reac2sidx
end

function build_active_set(pairs, L)
  isactive = zeros(Bool, L + 1)

  for pair in pairs
    reactant = pair[1]
    isactive[reactant] = true
  end

  return isactive
end

function build_lk2sidx_dict(compositions, reac2sidx, isactive, L)
  lk2sidx = OrderedDict{Tuple{Int,Int},Vector{Int}}()
  number_compositions = length(compositions)

  ### forming the l_kToSampleidx[] Dict now...

  for l in 1:(L + 1)

    # check if the species is active before doing anything
    if isactive[l]

      for k in 1:number_compositions

        ### get the composition associated with k, will need to change this to the appropriate dict call...
        comp = compositions[k]

        # get an empty vector to collect the sample classes a particle of type l w/ neighborhood k takes part in
        sample_classes = Int[]

        # iterate over every reactant pair [l, j]
        for j in 1:(L + 1)

          # get the number of neighbors of type j
          number_neighbors = comp[j]

          # check if there are any neighbors of the appropriate type
          if number_neighbors != 0

            # get the sample class, return 0 if [l, j] aren't reactant pairs
            s = get(reac2sidx, (l, j), 0)

            if s != 0

              # if [l, j] take part in a reaction, push the appropriate sample class to the vector
              push!(sample_classes, s + number_neighbors - 1)

            end
          end
        end

        # finally, check to see if [l, 0] is a reactant pair, meaning l takes part in an on-site reaction

        s = get(reac2sidx, (l, 0), 0)

        if s != 0
          # no need to correct s here
          push!(sample_classes, s)
        end

        # finally, put the sample_classes in the dictionary with key [l, k]
        lk2sidx[(l, k)] = sample_classes

      end
    end
  end
  return lk2sidx
end

##### dependency graphs #####

affect(rxn) = (rxn.input1, rxn.class)

# determine the types affected by a reaction
affected_types(rxn) = unique(setdiff([rxn.input1, rxn.input2, rxn.output1, rxn.output2], [0, 1]))

# enumerate the AffectedBy
function affected_by(rxn, sample_class, number_classes)
  # determine the types that can change
  affected_by_set = Tuple{Int,Int}[]

  for x in affected_types(rxn)
    for s in 1:number_classes
      for s_star in sample_class[(x, s)]
        affected_by_set = union(affected_by_set, [(x, s_star)])
      end
    end
  end

  return affected_by_set
end

count_pairwise(pairs) = count(x -> x[2] != 0, pairs)
count_onsite(pairs) = count(x -> x[2] == 0, pairs)

function rxnrxn_depgraph(::Iter, reaction::Vector{IPSReactionStruct}) where Iter <: DGIterStyle
  number_reactions = length(reaction)

  # determine the set of particle types
  ptypes = [ptype for rxn in reaction for ptype in (rxn.input1, rxn.input2, rxn.output1, rxn.output2)] |> unique

  # determine reactant pairs
  reactant_pairs = build_reactant_pairs(reaction)

  # count the number of unique particle types and neighborhood capacity
  L = setdiff(ptypes, [0, 1]) |> length

  ##### absolutely dirty hack
  nbmax = 0; rxn = reaction[nbmax+=1]

  search_tuple = (rxn.input1, rxn.input2, rxn.output1, rxn.output2)
  current_tuple = (rxn.input1, rxn.input2, rxn.output1, rxn.output2)

  while current_tuple == search_tuple
    rxn = reaction[nbmax+=1]
    current_tuple = (rxn.input1, rxn.input2, rxn.output1, rxn.output2)
  end

  # build the sample classes
  isactive = build_active_set(reactant_pairs, L)
  pair2class = build_pair_dictionary(reactant_pairs, nbmax)
  compositions = collect(Vector{Int64}, multiexponents(L + 1, nbmax))
  sample_class = build_lk2sidx_dict(compositions, pair2class, isactive, L)

  # count the number of sample classes
  number_classes = nbmax * count_pairwise(reactant_pairs) + count_onsite(reactant_pairs)

  # build the dependency graph
  deps_range = Vector{Tuple{Int,Int}}(undef, number_reactions)
  deps_value = Int[]

  idx = 1
  for j in eachindex(reaction)
    affected_by_j = affected_by(reaction[j], sample_class, number_classes)
    start = idx

    for k in eachindex(reaction)
      affect_k = affect(reaction[k])

      if affect_k in affected_by_j
        push!(deps_value, k)
        idx += 1
      end
    end
    stop = idx - 1

    deps_range[j] = (start, stop)
  end

  return DGVector{Iter}(deps_range, deps_value)
end
