export Network

type Network
    id::ASCIIString
    spcs::Dict{ASCIIString,Species}
    rxns::Dict{ASCIIString,Reaction}
    param::Dict{ASCIIString, Float64}
end

function Network()
  Network("untitled", Dict{ASCIIString,Species}(), Dict{ASCIIString,Reaction}(), Dict{ASCIIString,Float64}())
end

function Network(id)
  network = Network()
  set_id!(network, id)
  return network
end

function set_id!(network, id)
  network.id = id
  return network
end

function get_id(network)
  return network.id
end

# Species functions
function has_species(network, id)
  return haskey(network.spcs, id)
end

function has_species(network, s::Species)
  return has_species(network, s.id)
end

function add!(network, o::Species)
  id = o.id
  has_species(network, o) ? warn(o.id," already defined in ",network.id) : network.spcs[id] = o
  return network.spcs[id]
end

function delete!(network, o::Species)
  id = o.id
  has_species(network, o) ? delete!(network.spcs, id) : warn(o.id," not defined in ",network.id)
  return o
end

# function edit!(network, o::Species, values...)
#   if has_species(network, o)
#     d = Dict{Symbol, Any}()
#     for (key, value) in args
#       d[key] = value
#     end
#
#     s = network.spcs[o.id]
#     if haskey(d,:id); s.id = d[:id]::ASCIIString; end
#     if haskey(d,:pop); s.pop = d[:pop]::Int; end
#     if haskey(d,:istracked); s.istracked = d[:istracked]::Bool; end
#   else
#     error("Cannot edit ",s.id," in "network.id,": Undefined Species.")
#   end
#   return network
# end

# Reaction functions
function has_reaction(network, id)
  return haskey(network.rxns, id)
end

function has_reaction(network, r::Reaction)
  return has_reaction(network, r.id)
end

function add!(network, o::Reaction)
  id = o.id
  has_reaction(network, o) ? warn(o.id," already defined in ",network.id) : network.rxns[id] = o
  return network.rxns[id]
end

function delete!(network, o::Reaction)
  id = o.id
  has_reaction(network, o) ? delete!(network.rxns, id) : warn(o.id," not defined in ",network.id)
  return o
end
