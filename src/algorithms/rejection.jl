mutable struct RejectionMethod{S<:RatesCacheTrait} <: DecoupledSearchAlgorithm
  abstract_state::Vector{Tuple{Float64,Float64}}
  abstract_rates::Vector{Tuple{Float64,Float64}}
end

RatesCache(::RejectionMethod{S}) where S = S()

# TODO
# develop structs for intervals
# rewrite the search functions using the iteration interface

# OUTLINE

@inline function generate_jump(algorithm::RejectionMethod)
  j, u = accept_reject(RatesCache(algorithm), algorithm.abstract_rates, algorithm.total_rate)
  s = -log(u) / max(total_rate)
  
  return (j, s)
end

# rates should be an iterator over intervals_rates
@inline function accept_reject(ctype::S, abstract_rates, total_rate) where S
  while rejected
    # search jumps based on upper bounds
    j = search_jumps(ctype, rates, upper_bound(abstract_rates))

    rejected = rand() > lower_bound(abstract_rates, j) / upper_bound(abstract_rates, j)

    if rejected
      # evaluate propensity; how to get state and model here?
      rejected = rand() > rate_j / upper_bound(abstract_rates, j)
    end

    u *= rand()
  end

  return j, u
end

@inline function update!(algorithm::DirectMethod, state, model, j) = update_jump_rates!(algorithm, state, model)
  # check for states that have left interval
  # define a new interval for each vagabond state
  # check which propensity bounds must be updated based on a species-reaction dependency graph
end