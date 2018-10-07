mutable struct SortingDirectMethod{S<:RatesCacheTrait} <: DecoupledSearchAlgorithm
  rates::Vector{Float64}
  total_rate::Float64
  search_order::Vector{Int}
  search_index::Int
end

RatesCache(::SortingDirectMethod{S}) where S = S()

@inline function initialize!(algorithm::SortingDirectMethod, state, model)
  search_order = algorithm.search_order
  
  for i in eachindex(search_order)
      @inbounds search_order[i] = i
  end
  
  algorithm.search_index = 1
  
  update_jump_rates!(algorithm, state, model)

  return nothing
end

@inline function search_jumps(algorithm::SortingDirectMethod)
  idx = search_jump_rates(RatesCache(algorithm), algorithm.rates, algorithm.total_rate, algorithm.search_order)
  
  # record the index into the search order vector
  algorithm.search_index = idx
  
  # return the reaction index
  return @inbounds algorithm.search_order[idx]

  return nothing
end

@inline function update!(algorithm::SortingDirectMethod, state, model, j)
  # update dependent propensities
  update_jump_rates!(algorithm, state, model, j)
  
  search_order = algorithm.search_order
  search_index = algorithm.search_index
  
  # perform a sorting step
  if search_index > 1
      @inbounds tmp = search_order[search_index]
      
      @inbounds search_order[search_index]     = search_order[search_index - 1]
      @inbounds search_order[search_index - 1] = tmp
  end

  return nothing
end
