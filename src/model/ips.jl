struct IPSReactionStruct
  pairwise :: Bool
  input1   :: Int # center
  input2   :: Int # adjacent
  output1  :: Int # center
  output2  :: Int # adjacent
  class    :: Int
  rate     :: Float64
  klaw     :: KineticLaw
  # formula    :: Expr
  # parameter  :: Symbol
  # paramidx   :: Int

  function IPSReactionStruct(pairwise, i1, i2, o1, o2, class, rate, klaw)
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

    return new(pairwise, i1, i2, o1, o2, class, rate, klaw)
  end
end

function Base.show(io::IO, x::IPSReactionStruct)
  i1, i2 = x.input1, x.input2
  o1, o2 = x.output1, x.output2

  if x.pairwise
    str = "$(i1) + $(i2) --> $(o1) + $(o2)"
  else
    str = "$(i1) --> $(o1)"
  end

  print(io, str, "  class = ", x.class)
end

function Base.show(io::IO, m::MIME"text/plain", x::IPSReactionStruct)
  show(io, x)
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

# struct InteractingParticleSystem{DG,enumType}
struct InteractingParticleSystem{enumType}
  reactions::Vector{IPSReactionStruct}
  # rxn_rates::Vector{T}
  # dep_graph::DG
  # spc_graph::DG
  # rxn_graph::DG
  isactive::Vector{Bool}
  enumeration::enumType
end

function InteractingParticleSystem(reactions::Vector{IPSReactionStruct}, isactive, enumeration)
  num_reactions = length(reactions)

  # dep_graph = rxnrxn_depgraph(DGLazy(), reactions, isactive, enumeration)

  # return InteractingParticleSystem(reactions, rxn_rates, dep_graph, spc_graph, rxn_graph)
  # return InteractingParticleSystem(reactions, dep_graph, isactive, enumeration)
  return InteractingParticleSystem(reactions, isactive, enumeration)
end

function Base.summary(io::IO, x::InteractingParticleSystem)
  types = length(x.isactive) - 1
  str = length(types) > 1 ? "types" : "type"
  print(io, "Interacting Particle System w/ $(types) ", str)
end

function Base.show(io::IO, x::InteractingParticleSystem)
  summary(io, x)
  println()
  show(io, x.reactions)
end

function Base.show(io::IO, m::MIME"text/plain", x::InteractingParticleSystem)
  summary(io, x)
  println()
  show(io, m, x.reactions)
end


##### execute_jump! #####

execute_jump!(state, model::InteractingParticleSystem, j) = execute_jump!(state, model.reactions[j], model.enumeration)

function execute_jump!(lattice, reaction, enumeration)
  # xclass = lattice.xclass
  # xdiff  = lattice.xdiff

  # save counts in each class
  # for s in eachindex(xclass)
  #   xdiff[s] = false
  #   xclass[s] = length(lattice.classes[s].members)
  # end

  if ispairwise(reaction)
    # println("pairwise")
    execute_pairwise!(lattice, reaction, enumeration)
  else
    # println("onsite")
    execute_onsite!(lattice, reaction, enumeration)
  end

  # check which counts have changed
  # for s in eachindex(xclass)
  #   if xclass[s] != length(lattice.classes[s].members)
  #     xdiff[s] = true
  #   end
  # end

  return nothing
end

function execute_pairwise!(lattice, reaction, enumeration)
  # unpack information
  i1, i2 = get_reactant_types(reaction)
  j1, j2 = get_product_types(reaction)

  s = reaction.class
  # N = lattice.N

  # sample a particle of type i1 with nbclass k; i.e. sample class s
  x = sample_from_class(lattice, enumeration, s)
  #
  # if get_ptype(x) == 1
  #   println("execute reaction:")
  #   @show reaction
  #   println("center particle:")
  #   display(x)
  #   println("class members:")
  #   for z in lattice.site[enumeration.class[s]]
  #     print("  ")
  #     display(z)
  #   end
  #   error("sampled an open site")
  # end

  # sample a particle of type i2 from x's neighborhood
  y = sample_neighbor(x, lattice, enumeration, i2)
  # we use up the particles (k1, i1) and (k2, i2)

  # change the types for x and y
  i1 != j1 && transform!(x, j1)
  i2 != j2 && transform!(y, j2)

  # update neighborhood and sample classes for x and y
  update_classes_particle!(enumeration, x, i1, j1, i2, j2)
  update_classes_particle!(enumeration, y, i2, j2, i1, j1)

  # update the remaining particles surrounding x and y
  i1 != j1 && update_classes_neighbors!(lattice, enumeration, x, y, i1, j1)
  i2 != j2 && update_classes_neighbors!(lattice, enumeration, y, x, i2, j2)

  # we produced new particles so update the counts
  # N[get_neighbor_class(x), get_ptype(x)] += 1
  # N[get_neighbor_class(y), get_ptype(y)] += 1

  return nothing
end

function execute_onsite!(lattice, reaction, enumeration)
  # unpack reaction information
  i, _ = get_reactant_types(reaction)
  j, _ = get_product_types(reaction)

  s = reaction.class
  # N = lattice.N

  # sample a particle of type i with nbclass k; i.e. sample class s
  x = sample_from_class(lattice, enumeration, s)

  # we use up the x particle
  # N[get_neighbor_class(x), get_ptype(x)] -= 1

  # change the type for x
  transform!(x, j)

  # update neighborhood and sample classes for x and y
  update_classes_particle!(enumeration, x, i, j, 0, 0)

  # update the remaining particles surrounding x and y
  update_classes_neighbors!(lattice, enumeration, x, x, i, j)

  # update the counts
  # N[get_neighbor_class(x), get_ptype(x)] += 1

  return nothing
end

# dispatch to select reaction
function rate(model::InteractingParticleSystem, lattice, j)
  rate(model.reactions[j], lattice, model.enumeration)
end

# dispatch on rate law; reaction.klaw introduces type instability
function rate(reaction::IPSReactionStruct, lattice, enumeration)
  rate(reaction.klaw, reaction, lattice, enumeration)
end

# default implementation for rate
function rate(::DefaultIPSLaw, reaction::IPSReactionStruct, lattice, enumeration)
  # grab reactant indices
  class = enumeration.class
  s = reaction.class

  # grab per particle rate
  per_particle_rate = reaction.rate

  # grab the count for the center particle
  N = length(class[s]) # need to fix this

  return N * per_particle_rate
end
