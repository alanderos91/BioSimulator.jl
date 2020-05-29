##### KineticLaw / MassActionOrder* #####

abstract type KineticLaw end

struct MassAction  <: KineticLaw end
struct MichaelisMenten <: KineticLaw end

##### ReactionStruct #####

struct ReactionStruct{MA <: KineticLaw}
    reactants :: Vector{Tuple{Int,Int}}
    net_change :: Vector{Tuple{Int,Int}}
    paramidx :: Vector{Int}

    function ReactionStruct(law::MA, reactants, net_change, paramidx::Vector) where MA <: KineticLaw
        num_reactants = length(reactants)
        order = num_reactants > 0 ? sum(c for (_, c) in reactants) : 0

        !is_compatible_law(law, order, num_reactants) && throw(ArgumentError("reaction is not compatible with $(law) law"))

        return new{MA}(reactants, net_change, paramidx)
    end
end

# for test compatibility
ReactionStruct(law, reactants, net_change, paramidx::Real) = ReactionStruct(law, reactants, net_change, [paramidx])

function Base.show(io::IO, r::ReactionStruct{law}) where law <: KineticLaw
    # formula = r.formula
    # lhs = formula.args[1]
    # rhs = formula.args[2]
    #
    # print(io, string("  ", r.name, ": ", lhs, " -> ", rhs))
    print(io, "ReactionStruct")
end

is_compatible_law(::MassAction, order, num_reactants) = true
is_compatible_law(::MichaelisMenten, order, num_reactants) = (order ≤ 1)

function execute_jump!(x, r::ReactionStruct)
    net_change = r.net_change

    for v in net_change
        k, v_k = v
        @inbounds @fastmath x[k] += v_k
    end
end

##### rate functions for stochastic mass action kinetics #####

@inline @inbounds function rate(r::ReactionStruct{MassAction}, x, p)
    i = r.paramidx[1]
    total_rate = p[i]   # accumulates terms of the form (x_k - (j-1))
    prefactor = one(total_rate) # accumulates constants 1 / m_k!

    for (k, m) in r.reactants
        for j in 1:m
            total_rate *= (x[k] - (j-1))
            prefactor *= j
        end
    end

    total_rate /= prefactor

    return total_rate
end

@inline @inbounds function rate_derivative(r::ReactionStruct{MassAction}, x, p, i)
    total_rate = p[r.paramidx[1]]  # accumulates terms of the form (x_k - (j-1))
    prefactor = one(total_rate) # accumulates constants 1 / m_k!
    deriv_term = zero(total_rate)

    for (k, m) in r.reactants
        for j in 1:m
            total_rate *= (x[k] - (j-1))
            prefactor *= j

            if k == i
                deriv_term += 1 / (x[k] - (j-1))
            end
        end
    end

    total_rate = total_rate / prefactor * deriv_term

    return total_rate
end

##### rate functions for stochastic Michaelis-Menten kinetics #####
@inline @inbounds function rate(r::ReactionStruct{MichaelisMenten}, x, p)
    i1 = r.paramidx[1]
    i2 = r.paramidx[2]

    V = p[i1]
    K = p[i2]
    k = r.reactants[1][1]

    total_rate = V * x[k] / (K + x[k])

    return total_rate
end

##### ReactionSystem #####
ReactionLike = Union{ReactionStruct{MassAction},ReactionStruct{MichaelisMenten}}

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
    rxn_rates = Float64[]

    build_reactions!(reactions, rxn_rates, model)

    dep_graph = rxnrxn_depgraph(DGLazy(), model)
    spc_graph = spcrxn_depgraph(DGLazy(), model)
    rxn_graph = rxnspc_depgraph(DGLazy(), model)

    return ReactionSystem(reactions, rxn_rates, dep_graph, spc_graph, rxn_graph)
end

function Base.summary(io::IO, r::ReactionSystem)
    print(io, "Well-Mixed Reaction System")
end

function Base.show(io::IO, r::ReactionSystem)
    summary(io, r)
    for (j, reaction) in enumerate(r.reactions)
        print(io, "\n", reaction)
    end
end

##### convenience functions #####

@inline function execute_jump!(x, rxn::ReactionSystem, j)
    @inbounds execute_jump!(x, rxn.reactions[j])
end

@inline @inbounds rate(rxn::ReactionSystem, x, j) = rate(rxn.reactions[j], x, rxn.rxn_rates)

@inline @inbounds function rate_derivative(rxn, x, i, j)
    rate_derivative(rxn.reactions[j], x, rxn.rxn_rates, i)
end

function netstoichiometry(rxn::ReactionSystem, num_species, num_reactions)
    V = zeros(Int, num_species, num_reactions)
    reactions = rxn.reactions
    for j in eachindex(reactions)
        r = reactions[j]
        for (k, v) in r.net_change
            V[k, j] = v
        end
    end
    return V
end

##### helper functions for building a ReactionSystem #####

function build_reactions!(rxn_set, rxn_rates, model)
    species   = species_list(model)
    reactions = reaction_list(model)

    indexmap = OrderedDict(key => i for (i, key) in enumerate(keys(species)))

    pidx = 1
    for (j, r) in enumerate(values(reactions))
        name = r.identifier
        formula = r.origex
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

        if r.law == :mass_action
            klaw = MassAction()
            paramidx = [pidx]
            pidx += 1
            push!(rxn_rates, rxn_rate[1])
        elseif r.law == :michaelis_menten
            klaw = MichaelisMenten()
            paramidx = [pidx, pidx+1]
            pidx += 2
            push!(rxn_rates, rxn_rate[1])
            push!(rxn_rates, rxn_rate[2])
        else
            error("unknown kinetic law $(r.law)")
        end

        rxn_set[j] = ReactionStruct(klaw, rtuples, net_change, paramidx)
    end

    return rxn_set
end

function get_kinetic_law(rtuples)
    # num_reactants = length(rtuples)
    # order = num_reactants > 0 ? sum(c for (_, c) in rtuples) : 0
    return MassAction()
end

##### extras #####

function execute_leap!(state, stoichiometry, number_jumps)
    state += stoichiometry * number_jumps
end
