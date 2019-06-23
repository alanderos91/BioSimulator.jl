struct SampleClassEnumeration
  class::Vector{Int}
  composition::Vector{Vector{Int}}
  reactant_to_class::OrderedDict{Tuple{Int,Int},Int}
  pair_to_classes::OrderedDict{Tuple{Int,Int},Vector{Int}}
  dummy_composition::Vector{Int}
end

function SampleClassEnumeration(composition, reactant_to_class, pair_to_classes)
  class = Int[]
  dummy_composition = zero(composition[1])

  SampleClassEnumeration(
    class,
    composition,
    reactant_to_class,
    pair_to_classes,
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
