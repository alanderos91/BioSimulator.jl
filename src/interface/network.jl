import Base.show

type Network
  id::ASCIIString
  species::Dict{Symbol,Species}
  reactions::Dict{Symbol,Reaction}
  parameters::Dict{Symbol,Parameter}

  function Network(id)
    s = Dict{Symbol,Species}()
    r = Dict{Symbol,Reaction}()
    p = Dict{Symbol,Parameter}()
    return new(id, s, r, p)
  end
end

function Base.show(io::IO, x::Network)
  @printf io "[ Model: %s ]\n" x.id
  @printf io " no. species:    %d\n" length(x.species)
  @printf io " no. reactions:  %d\n" length(x.reactions)
  @printf io " no. parameters: %d"   length(x.parameters)
end

function add_object!(model, object, fieldname)
  dict = getfield(model, fieldname)
  id   = object.id
  setindex!(dict, object, id)
  return dict
end

function (<=)(model::Network, object::Species)
  add_object!(model, object, :species)
end

function (<=)(model::Network, object::Reaction)
  reactants = object.reactants
  products  = object.products
  try
    validate(reactants, model)
    validate(products,  model)
    add_object!(model, object, :reactions)
  catch ex
    rethrow(ex)
  end
end

function (<=)(model::Network, object::Parameter)
  add_object!(model, object, :parameters)
end

function validate(participants, model)
  species_dict = model.species
  for species in keys(participants)
    if !haskey(species_dict, species)
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

function (>=)(model::Network, x::Pair{Symbol,Symbol})
  fieldname = x.first
  key = x.second

  if fieldname == :species
    # remove all species refs in reactions?
  end
  
  rmv_object!(model, key, fieldname)
end
