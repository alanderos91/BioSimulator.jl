immutable Model{T<:AbstractReactionSystem}
    id::ASCIIString
    Xt::Vector{Int}
    rs::T

    parameters::Dict{Symbol,Parameter}

    X₀::Vector{Int}
end

##### accessors #####
id(m::Model)         = m.id
species(m::Model)    = m.Xt
reactions(m::Model)  = m.rs
parameters(m::Model) = m.parameters
initial(m::Model)    = m.X₀

function reset!(m::Model)
  copy!(m.Xt, m.X₀)
  m.rs.a.a0 = Inf
  #fill!(m.rv, 0.0)
end

function make_species_arr(dict::Dict{Symbol,Species})
    Xt     = Array(Int,    length(dict))
    id     = Array(Symbol, length(dict))
    id2ind = Dict{Symbol,Int}()

    i = 1
    for (key, s) in dict
        Xt[i] = s.population
        id[i] = s.id
        id2ind[key] = i
        i = i + 1
    end
    return Xt, id, id2ind
end
