mutable struct DirectMethod{S<:RatesCacheTrait} <: DecoupledSearchAlgorithm
  rates::Vector{Float64}
  total_rate::Float64
end

RatesCache(::DirectMethod{S}) where S = S()

@inline search_jumps(algorithm::DirectMethod) = search_jump_rates(RatesCache(algorithm), algorithm.rates, algorithm.total_rate)

@inline update!(algorithm::DirectMethod, state, model, j) = update_jump_rates!(algorithm, state, model)
