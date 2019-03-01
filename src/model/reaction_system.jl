##### KineticLaw / MassActionOrder* #####

abstract type KineticLaw end

struct MassActionOrder0  <: KineticLaw end
struct MassActionOrder1  <: KineticLaw end
struct MassActionOrder2A <: KineticLaw end
struct MassActionOrder2B <: KineticLaw end

##### ReactionStruct #####

struct ReactionStruct{MA <: KineticLaw}
  reactants :: Vector{Tuple{Int,Int}}
  net_change :: Vector{Tuple{Int,Int}}
  paramidx :: Int

  function ReactionStruct(law::MA, reactants, net_change, paramidx) where MA <: KineticLaw
    num_reactants = length(reactants)
    order = num_reactants > 0 ? sum(c for (_, c) in reactants) : 0

    !is_compatible_law(law, order, num_reactants) && throw(ArgumentError("reaction is not compatible with $(law) law"))

    return new{MA}(reactants, net_change, paramidx)
  end
end

is_compatible_law(::MassActionOrder0,  order, num_reactants) = order == 0
is_compatible_law(::MassActionOrder1,  order, num_reactants) = order == 1
is_compatible_law(::MassActionOrder2A, order, num_reactants) = order == 2 && num_reactants == 2
is_compatible_law(::MassActionOrder2B, order, num_reactants) = order == 2 && num_reactants == 1

function execute_jump!(x, r::ReactionStruct)
    net_change = r.net_change
    
    for v in net_change
        k, v_k = v
        @inbounds @fastmath x[k] += v_k
    end
end

##### rate functions for stochastic mass action kinetics #####

@inline @inbounds function rate(r::ReactionStruct{MassActionOrder0}, x, p)
  i = r.paramidx
  
  return p[i]
end

@inline @inbounds function rate(r::ReactionStruct{MassActionOrder1}, x, p)
  i = r.paramidx
  k, _ = r.reactants[1]
  
  return p[i] * x[k]
end

@inline @inbounds function rate(r::ReactionStruct{MassActionOrder2A}, x, p)
  i = r.paramidx
  k1, _ = r.reactants[1]
  k2, _ = r.reactants[2]
  
  return p[i] * x[k1] * x[k2]
end

@inline @inbounds function rate(r::ReactionStruct{MassActionOrder2B}, x, p)
  i = r.paramidx
  k, _ = r.reactants[1]
  
  return 1//2 * x[k] * (x[k] - 1) * p[i]
end

##### type union for heterogeneous ReactionStruct arrays #####

ReactionLike = Union{
    ReactionStruct{MassActionOrder0},
    ReactionStruct{MassActionOrder1},
    ReactionStruct{MassActionOrder2A},
    ReactionStruct{MassActionOrder2B}
  }

##### ReactionSystem #####

struct ReactionSystem{R,DG<:DependencyGraph}
  reactions::Vector{ReactionLike}
  rxn_rates::R
  dep_graph::DG
  spc_graph::DG
  rxn_graph::DG
end

function ReactionSystem(model::Network)
  num_reactions = number_reactions(model)
  
  reactions = Vector{ReactionLike}(undef, num_reactions)
  rxn_rates = zeros(num_reactions)

  build_reactions!(reactions, rxn_rates, model)

  dep_graph = rxnrxn_depgraph(DGLazy(), model)
  spc_graph = spcrxn_depgraph(DGLazy(), model)
  rxn_graph = rxnspc_depgraph(DGLazy(), model)

  return ReactionSystem(reactions, rxn_rates, dep_graph, spc_graph, rxn_graph)
end

##### convenience functions #####

@inline function execute_jump!(x, rxn::ReactionSystem, j)
    @inbounds execute_jump!(x, rxn.reactions[j])
end

@inline @inbounds rate(rxn::ReactionSystem, x, j) = rate(rxn.reactions[j], x, rxn.rxn_rates)

##### helper functions for building a ReactionSystem #####

function build_reactions!(rxn_set, rxn_rates, model)
  species   = species_list(model)
  reactions = reaction_list(model)

  indexmap = OrderedDict(key => i for (i, key) in enumerate(keys(species)))

  for (j, r) in enumerate(values(reactions))
    reactants = r.reactants
    products  = r.products
    rxn_rate  = r.rate

    rtuples = Tuple{Int,Int}[(indexmap[s], c) for (s, c) in reactants]
    net_change = Tuple{Int,Int}[]

    for s in keys(species)
      change = get(products, s, 0) - get(reactants, s, 0)
      if change != 0
        push!(net_change, (indexmap[s], change))
      end
    end
    L = get_kinetic_law(rtuples)
        
    rxn_set[j] = ReactionStruct(L(), rtuples, net_change, j)
    rxn_rates[j] = rxn_rate
  end

  return rxn_set
end

function get_kinetic_law(rtuples)
    num_reactants = length(rtuples)
    order = num_reactants > 0 ? sum(c for (_, c) in rtuples) : 0
    
    if order == 0
        MassActionOrder0
    elseif order == 1
        MassActionOrder1
    elseif order == 2 && num_reactants == 2
        MassActionOrder2A
    elseif order == 2 && num_reactants == 1
        MassActionOrder2B
    else
        error("reaction is not elementary or mass action")
    end
end
