mutable struct FirstReactionMethod <: CoupledSearchAlgorithm
  rates::Vector{Float64}
  total_rate::Float64
end

RatesCache(::FirstReactionMethod) = HasRates()

@inline search_jumps(algorithm::FirstReactionMethod) = search_jump_times(RatesCache(algorithm), algorithm.rates)

@inline update!(algorithm::FirstReactionMethod, state, model, j) = update_jump_rates!(algorithm, state, model)
