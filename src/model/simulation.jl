export Simulation

immutable Simulation
    id::ASCIIString
    initial::Vector{Species}

    state::Vector{Species}
    rxns::Vector{Reaction}
    param::Dict{ASCIIString, Float64}
end

function Simulation(x::Network)
    return Simulation(x.id, deepcopy(x.spcs), deepcopy(x.spcs), deepcopy(x.rxns), deepcopy(x.param))
end
