function compute_net_stoichiometry!(V, U, reactions, id2ind)
  j = 1
  for (key, r) in reactions
    reactants = r.reactants
    products  = r.products

    for (id, i) in id2ind
      x = get(reactants, id, 0)
      y = get(products,  id, 0)
      V[i, j] = y - x
      U[i, j] = x
    end
    j = j + 1
  end

  return V, U
end

function compute_scaled_rates!(k, reactions)
  j = 1
  for (key, reaction) in reactions
    rate = reaction.rate
    temp = 1.0
    reactants = reaction.reactants

    for (reactant, coefficient) in reactants
      temp = temp * coefficient
    end

    k[j] = rate / temp

    j = j + 1
  end

  return k
end

function make_dependency_graph(reactions)
  dg = Vector{Int}[ [] for i in 1:length(reactions) ]
  make_dependency_graph!(dg, reactions)

  return dg
end

function make_dependency_graph!(dg, reactions_dict)
  reactions = values(reactions_dict)

  for (i, reaction) in enumerate(reactions)
    for (j, other) in enumerate(reactions)
      affected_species = union(keys(reaction.reactants), keys(reaction.products))
      other_reactants  = keys(other.reactants)

      if !isempty(intersect(affected_species, other_reactants))
        push!(dg[i], j)
      end
    end
    #dg[i] = unique(dg[i]) do we need this?
  end

  return dg
end

function make_species_vector(dict::Associative{Symbol,Species})
  X0     = Array{Int}(length(dict))
  id     = Array{Symbol}(length(dict))
  id2ind = Dict{Symbol,Int}()

  i = 1
  for (key, s) in dict
    X0[i] = s.population
    id[i] = s.id
    id2ind[key] = i
    i = i + 1
  end

  return X0, id, id2ind
end
