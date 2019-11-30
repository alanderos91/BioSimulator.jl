mutable struct NextReactionMethod{PQ<:PriorityQueue} <: CoupledSearchAlgorithm
  rates::Vector{Float64}
  total_rate::Float64
  pq::PQ
end

RatesCache(::NextReactionMethod) = HasRates()

@inline function initialize!(algorithm::NextReactionMethod, state, model)
  total_rate = zero(algorithm.total_rate)
  
  for j in eachindex(algorithm.rates)
    # get rate
    rate_j = rate(model, state, j)
    
    # compute putative reaction time
    t_j = rate_j > 0 ? randexp() / rate_j : Inf
    
    # store rate and reaction time
    @inbounds algorithm.rates[j] = rate_j
    @inbounds algorithm.pq[j]    = t_j
    
    total_rate = total_rate + rate_j
  end
  
  algorithm.total_rate = total_rate
end

# NRM returns absolute times, not increments
get_new_time(::NextReactionMethod, t, s) = s

@inline function search_jumps(algorithm::NextReactionMethod)
  j, s = peektop(algorithm.pq)
  
  return (j, s)
end

@inline function update!(algorithm::NextReactionMethod, state, model, k)
  rates = algorithm.rates
  pq = algorithm.pq
  t  = pq[k]
  
  for j in dependents(model.dep_graph, k) # returns a range
    if j != k
      # compute new rate
      rate_j = rate(model, state, j)
      
      @inbounds t_j = pq[j]

      # check for infinite jump times
      if isfinite(t_j)
        @inbounds pq[j] = rate_j > 0 ? t + (rates[j] / rate_j) * (t_j - t) : Inf
      else
        @inbounds pq[j] = rate_j > 0 ? t + randexp() / rate_j : Inf
      end
      
      # update total rate
      @inbounds algorithm.total_rate += rate_j - rates[j]
      
      # update reaction rate
      @inbounds rates[j] = rate_j
    end
  end
  
  # update the current reaction's propensity and firing time
  rate_k = rate(model, state, k)
  @inbounds pq[k] = rate_k > 0 ? t + randexp() / rate_k : Inf
  algorithm.total_rate += rate_k - rates[k]
  @inbounds rates[k] = rate_k
  
  return algorithm.total_rate
end
