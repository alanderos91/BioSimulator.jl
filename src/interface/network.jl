import Base.show

"""
```
Network(id)
```

Constructs an empty `Network` that models a system of interacting particles. A particle is represented by a `Species`; an interaction is represented by a `Reaction`. Model parameters (e.g. an intrinsic birth rate) are represented by a `Parameter` type.

Adding a `Species`, `Reaction`, or `Parameter` is done using the `<=` operator:

```
m <= Species("X", 100) # Adds a Species
m <= Reaction("birth", 2.0, :(X --> X + X)) # Adds a Reaction
```

A fully specified `Network` model is simulated by calling `simulate`:

```
simulate(m, SSA(4.0))
```

A `Network` is visualized by calling `petrinet`:

```
petrinet(m)
```

See also: `Species`, `Reaction`, `Parameter`, `simulate`, `petrinet`

### Arguments
`id`: An string identifier for the `Network`.
"""
type Network
  id :: UTF8String

  species_list   :: Dict{Symbol,Species}
  reaction_list  :: Dict{Symbol,Reaction}

  function Network(id)
    s = Dict{Symbol,Species}()
    r = Dict{Symbol,Reaction}()

    return new(UTF8String(id), s, r)
  end
end

species_list(x::Network)   = x.species_list
reaction_list(x::Network)  = x.reaction_list

n_species(x::Network)    = length(species_list(x))
n_reactions(x::Network)  = length(reaction_list(x))

function Base.show(io::IO, x::Network)
  print(io, "[ Model: $(x.id) ]\n")
  print(io, " no. species:    $(n_species(x))\n")
  print(io, " no. reactions:  $(n_reactions(x))")
end

function add_object!(model, object, fieldname)
  dict = getfield(model, fieldname)
  id   = object.id
  setindex!(dict, object, symbol(id))
  return dict
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
    add_object!(model, object, :reaction_list)
  catch ex
    rethrow(ex)
  end
end

function validate(participants, model)
  species_dict = species_list(model)
  for species in keys(participants)
    if !haskey(species_dict, symbol(species))
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

# function (>=)(model::Network, x::Pair{Symbol,UTF8String})
#   fieldname = x.first
#   key = x.second
#
#   if fieldname == :species_list || :species
#     # remove all species refs in reactions?
#   end
#
#   rmv_object!(model, key, fieldname)
# end
