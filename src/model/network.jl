export Network

type Network
    id::ASCIIString
    spcs::Vector{Species}
    rxns::Vector{Reaction}
    param::Dict{ASCIIString, Float64}
end
