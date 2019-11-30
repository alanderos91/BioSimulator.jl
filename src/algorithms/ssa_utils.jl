abstract type RatesCacheTrait end

struct HasRates <: RatesCacheTrait end

# populate the rates vector with each jump's propensity
@inline function update_jump_rates!(::HasRates, rates, total_rate, state, model)
  total_rate = zero(total_rate)
  
  for j in eachindex(rates)
      rate_j = rate(model, state, j)
      @inbounds rates[j] = rate_j
      @fastmath total_rate += rate_j
  end
  
  return total_rate
end

# populate the rates vector with each reaction channel's propensity evaluated at x,
# given that reaction R_k occured
@inline function update_jump_rates!(::HasRates, rates, total_rate, state, model, k)
  for j in dependents(model.dep_graph, k)
      rate_j = rate(model, state, j)
      @inbounds @fastmath total_rate += rate_j - rates[j]
      @inbounds rates[j] = rate_j
  end
  
  return total_rate
end

@inline function update_jump_rates!(::HasRates, rates, total_rate, state, model::InteractingParticleSystem, unused)
  xdiff = state.xdiff
  reactions = model.reactions

  for s in eachindex(xdiff)
    if xdiff[s]
      for j in eachindex(rates)
        rxn = reactions[j]

        if s == rxn.class
          rate_j = rate(model, state, j)
          @inbounds @fastmath total_rate += rate_j - rates[j]
          @inbounds rates[j] = rate_j
        end
      end
    end
  end

  return total_rate
end

# linear search
@inline function search_jump_rates(::HasRates, rates, total_rate)
  u = rand() * total_rate
  c = zero(total_rate)
  j = 0

  for rate in rates
    c += rate
    j += 1
    c >= u && break
  end

  return j
end

struct HasSums <: RatesCacheTrait end

# populate the rates vector with partial sums of the reaction channels' propensities evaluated at x
@inline function update_jump_rates!(::HasSums, rates, total_rate, state, model)
    rates[1] = rate(model, state, 1)
    
    for j in 2:length(rates)
        @fastmath rates[j] = rates[j - 1] + rate(model, state, j)
    end
    
    return rates[end]
end

# binary search on jump rates
@inline function search_jump_rates(::HasSums, rates, total_rate)
    return searchsortedfirst(rates, rand() * total_rate)
end

# search on jump times
@inline function search_jump_times(::HasRates, rates)
    n = length(rates)
    
    j = 1
    @inbounds s = randexp() / rates[1]
    
    for k in 2:n
        @inbounds w = randexp() / rates[k]
        if w < s
            j = k
            s = w
        end
    end
    
    return (j, s)
end