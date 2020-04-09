"""
```
Network(id)
```

Construct an empty `Network` representing a system of interacting particles. A particle is represented by a `Species`, and an interaction is represented by a `Reaction`.

Add a `Species` or `Reaction` using `<=`:

```
m <= Species("X", 100)
m <= Reaction("birth", 2.0, "X --> X + X")
```
"""
struct Network
  identifier :: Symbol

  species_list  :: OrderedDict{Symbol,Species}
  reaction_list :: OrderedDict{Symbol,Reaction}
  dep_graph :: OrderedDict{Symbol,Vector{Symbol}}

  function Network(id)
    s  = OrderedDict{Symbol,Species}()
    r  = OrderedDict{Symbol,Reaction}()
    dep_graph = OrderedDict{Symbol,Vector{Symbol}}()

    return new(Symbol(id), s, r, dep_graph)
  end
end

"""
```
species_list(x::Network)
```

Retrieve `Species` from a `Network`.
"""
species_list(x::Network) = x.species_list

"""
```
reaction_list(x::Network)
```

Retrieve `Reaction`s from a `Network`.
"""
reaction_list(x::Network) = x.reaction_list

get_species(m::Network, key::AbstractString) = get_species(m, Symbol(key))

function get_species(m::Network, key::Symbol)
  s = get(m.species_list, key, nothing)

  if isnothing(s)
    error("""
    Species $(key) not defined.
    Add it with `m <= Species("$(key)")`.
    """)
  end

  return s
end

get_reaction(m::Network, key::AbstractString) = get_reaction(m, Symbol(key))

function get_reaction(m::Network, key::Symbol)
  r = get(m.reaction_list, key, nothing)

  if isnothing(r)
    error("""
    Reaction $(key) not defined.
    Add it with `m <= Reaction("$(key)", rate, formula)`.
    """)
  end

  return r
end

number_species(x::Network)   = length(species_list(x))
number_reactions(x::Network) = length(reaction_list(x))

function Base.show(io::IO, x::Network)
  print(io, "[ Model: $(x.identifier) ]\n")
  print(io, "  no. species:    $(number_species(x))\n")
  print(io, "  no. reactions:  $(number_reactions(x))")
end

function add_object!(model, object, fieldname)
  dict = getfield(model, fieldname)
  id   = object.identifier
  setindex!(dict, object, Symbol(id))
  return model
end

function (<=)(model::Network, object::Species)
  add_object!(model, object, :species_list)
end

function (<=)(model::Network, object::Reaction)
  reactants = object.reactants
  products  = object.products

  try
    validate(reactants, model)
    validate(products,  model)

    # add the reaction to the list of objects
    add_object!(model, object, :reaction_list)

    id = Symbol(object.identifier)
    dep_graph = model.dep_graph

    if haskey(dep_graph, id)
      # if the key already exists, empty the list of dependents
      empty!(dep_graph[id])

      # and remove the reaction as a dependency
      for (key, r) in reaction_list(model)
        list = dep_graph[key]
        ix = findfirst(list, id)
        ix > 0 && deleteat!(list, ix)
      end
    else
      # if the key does not exist, make a new vector for dependents
      dep_graph[id] = Symbol[]
    end

    # update the dependency graph
    update_dep_graph!(model, object)

  catch ex
    rethrow(ex)
  end
end

function validate(participants, model)
  species_dict = species_list(model)
  for species in keys(participants)
    if !haskey(species_dict, Symbol(species))
      error("$(species) is not defined.")
    end
  end
end

function rmv_object!(model, key, fieldname)
  dict = getfield(model, fieldname)
  if haskey(dict, key)
    delete!(dict, key)
  else
    error("$(key) not found in $(fieldname)")
  end
  return model
end

function update_dep_graph!(model, rj)
  id = Symbol(rj.identifier)

  rj_affects    = rj.affects
  rj_affectedby = rj.affectedby

  for (key, rk) in reaction_list(model)
    rk_affects    = rk.affects
    rk_affectedby = rk.affectedby


    # add edge R_j --> R_k
    if !isempty(intersect(rj_affectedby, rk_affects))
      push!(model.dep_graph[id], key)
    end

    # add edge R_k --> R_j
    if rj != rk && !isempty(intersect(rk_affectedby, rj_affects))
      push!(model.dep_graph[key], id)
    end
  end

  return model
end
