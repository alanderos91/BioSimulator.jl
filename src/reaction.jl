type Reaction
    id::ASCIIString
    rate::ASCIIString
    propensity::Float64
    pre::Vector{Int}
    post::Vector{Int}
end