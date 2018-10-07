mutable struct EnhancedDirectMethod{S<:RatesCacheTrait} <: DecoupledSearchAlgorithm
  rates::Vector{Float64}
  total_rate::Float64
end

RatesCache(::EnhancedDirectMethod{S}) where S = S()

@inline search_jumps(algorithm::EnhancedDirectMethod) = search_jump_rates(RatesCache(algorithm), algorithm.rates, algorithm.total_rate)

@inline update!(algorithm::EnhancedDirectMethod, state, model, j) = update_jump_rates!(algorithm, state, model, j)
