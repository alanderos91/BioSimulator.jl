import Base.show

type Network
  id::ASCIIString
  species::Dict{ASCIIString,Species}
  reactions::Dict{ASCIIString,ReactionDef}
  parameters::Dict{ASCIIString,Float64}

  function Network(id)
    s = Dict{ASCIIString,Species}()
    r = Dict{ASCIIString,ReactionDef}()
    p = Dict{ASCIIString,Float64}()
    return new(id, s, r, p)
  end
end

function Base.show(io::IO, x::Network)
  print(io, "[ Model: ",x.id," ]\n")
  print(io, " Species:    ", length(x.species),"\n")
  print(io, " Reactions:  ", length(x.reactions),"\n")
  print(io, " Parameters: ", length(x.parameters))
end

function add_species!(network, sname, initial::Integer=0, istracked::Bool=true)
  d = getfield(network, :species)
  if haskey(d, sname)
    warn(string("Species ",sname," already defined."))
  else
    s = Species(sname, initial, istracked)
    setindex!(d, s, sname)
  end
  return d
end

function rmv_species!(network, sname)
  d = getfield(network, :species)
  if haskey(d, sname)
    delete!(d, sname)
    _remove_all_instances!(network, sname)
  else
    warn(string("Species ",sname," not defined."))
  end
  return d
end

function add_reaction!(network, rname, pname)
  d = getfield(network, :reactions)
  if haskey(d, rname)
    warn(string("Reaction ",rname," already defined."))
  else
    r = ReactionDef(rname, pname)
    setindex!(d, r, rname); add_parameter!(network, pname)
  end
  return d
end

function rmv_reaction!(network, rname)
  d = getfield(network, :reactions)
  if haskey(d, rname)
    r = getindex(d, rname)
    delete!(d, rname); rmv_parameter!(network, r.rate)
  else
    warn(string("Reaction ",rname," not defined."))
  end
  return d
end

function add_reactant!(network, rname, sname, coeff::Integer=1)
  sd = getfield(network, :species)
  rd = getfield(network, :reactions)
  if haskey(sd, sname)
    if haskey(rd, rname)
      r = getindex(rd, rname)
      add_participant!(r, sname, coeff, :reactants)
    else
      warn(string("Reaction ",rname," not defined."))
    end
  else
    warn(string("Species ",sname," not defined."))
  end
  return rd
end

function rmv_reactant!(network, rname, sname)
  rd = getfield(network, :reactions)
  if haskey(rd, rname)
    r = getindex(rd, rname)
    rmv_participant!(r, sname, :reactants)
  else
    warn(string("Reaction ",rname," not defined."))
  end
  return rd
end

function add_product!(network, rname, sname, coeff::Integer=1)
  sd = getfield(network, :species)
  rd = getfield(network, :reactions)
  if haskey(sd, sname)
    if haskey(rd, rname)
      r = getindex(rd, rname)
      add_participant!(r, sname, coeff, :products)
    else
      warn(string("Reaction ",rname," not defined."))
    end
  else
    warn(string("Species ",sname," not defined."))
  end
  return rd
end

function rmv_product!(network, rname, sname)
  rd = getfield(network, :reactions)
  if haskey(rd, rname)
    r = getindex(rd, rname)
    rmv_participant!(r, sname, :products)
  else
    warn(string("Reaction ",rname," not defined."))
  end
  return rd
end

function add_parameter!(network, pname, value::Float64=0.0)
  d = getfield(network, :parameters)
  if !haskey(d, pname)
    setindex!(d, value, pname)
  end
end

function set_parameter!(network, pname, value)
  d = getfield(network, :parameters)
  setindex!(d, value, pname)
end

function rmv_parameter!(network, pname)
  d = getfield(network, :parameters)
  delete!(d, pname)
end

function _remove_all_instances!(network, sname)
  rd = getfield(network, :reactions)
  for r in values(rd)
    rmv_participant!(r, sname, :reactants)
    rmv_participant!(r, sname, :products)
  end
end
