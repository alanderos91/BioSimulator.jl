import Base.show

type Network
  id::ASCIIString
  species::Dict{Symbol,Species}
  reactions::Dict{Symbol,Reaction}
  parameters::Dict{Symbol,Parameter}

  function Network(id)
    s = Dict{ASCIIString,Species}()
    r = Dict{ASCIIString,Reaction}()
    p = Dict{ASCIIString,Parameter}()
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
  if haskey(dict, id)
    warn()
  else
    setindex!(dict, object, id)
  end
  return dict
end

function (<=)(model::Network, object::Species)
  add_object!(model, object, :species)
end

function (<=)(model::Network, object::Reaction)
  # TODO: Check if species list is empty
  # TODO: Check if species in reaction exist in model
  add_object!(model, object, :reactions)
end

function (<=)(model::Network, object::Parameter)
  add_object!(model, object, :parameters)
end

# function rmv_object!(model, object, fieldname)
#   dict = getfield(model, fieldname)
#   if haskey(dict, id)
#     delete!(dict, id)
#   else
#     warn()
#   end
#   return dict
# end
#
# function (>=)(model::Network, object::Species)
#   # Remove references to species
#   rmv_object!(model, object, :species)
# end
#
# function (>=)(model::Network, object::Reaction)
#   rmv_object!(model, object, :reaction)
# end
#
# function (>=)(model::Network, object::Parameter)
#   rmv_object!(model, object, :species)
# end
