immutable Model
    id::Symbol
    Xt::Vector{Int}
    rv::ReactionVector

    parameters::Dict{Symbol,Parameter}

    X₀::Vector{Int}
end

##### accessors #####
id(m::Model)         = m.id
species(m::Model)    = m.Xt
reactions(m::Model)  = m.rv
parameters(m::Model) = m.parameters
initial(m::Model)    = m.X₀

function reset!(m::Model)
  copy!(m.Xt, m.X₀)
  m.rv.intensity = Inf
  #fill!(m.rv, 0.0)
end

compute_propensities!(m) = compute_intensities!(reactions(m), species(m), parameters(m))

function fire_reaction!(m::Model, r::ReactionChannel)
    Xt = species(m)
    @inbounds for i in eachindex(Xt)
        Xt[i] = Xt[i] + increment(r, i)
    end
end

function fire_reactions!(m, events)
    rxns = reactions(m)
    @inbounds for i in eachindex(rxns)
        fire_reaction!(m, reaction(rxns, i), events[i])
    end
end

function fire_reaction!(m, rxn, k)
    Xt = species(m)
    @inbounds for i in eachindex(Xt)
        Xt[i] = Xt[i] + k * increment(rxn, i)
    end
end

function make_species_arr(dict::Dict{Symbol,Species})
    Xt     = Array(Int,    length(dict))
    id     = Array(Symbol, length(dict))
    id2ind = Dict{Symbol,Int}()

    i = 1
    for (key, s) in dict
        Xt[i]       = s.population
        id[i]       = s.id
        id2ind[key] = i
        i = i + 1
    end
    return Xt, id, id2ind
end

function make_reaction_arr(dict::Dict{Symbol,Reaction}, id2ind)
    rxn   = Array(ReactionChannel, length(dict))
    id    = Array(Symbol,          length(dict))

    j = 1
    for (key, r) in dict
        rxn[j]     = ReactionChannel(r, id2ind)
        id[j]      = r.id
        j = j + 1
    end
    return ReactionVector(rxn), id
end
