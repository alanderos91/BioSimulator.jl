mutable struct RejectionMethod{S<:RatesCacheTrait} <: DecoupledSearchAlgorithm

end

RatesCache(::DirectMethod{S}) where S = S()

# TODO
# develop structs for intervals

# OUTLINE

@inline function generate_jump(algorithm::RejectionMethod)
  j, u = accept_reject(...)
  s = -log(u) / max(total_rate)
  
  return (j, s)
end

@inline function accept_reject()
  while rejected
    # search jumps based on upper bounds
    j = search_jumps(...)

    rejected = rand() > min(rates[j]) / max(rates[j])

    if rejected
      # evaluate propensity
      rejected = rand() > rates[j] / max(rates[j])
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