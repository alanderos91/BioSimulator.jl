mutable struct PoissonLeapAlgorithm{LF,C} <: UnsafeLeapAlgorithm
  rates::Vector{Float64}
  total_rate::Float64
  leap_formula::LF
  cache::C
end

function generate_leap!(v, algorithm::PoissonLeapAlgorithm)
  s = tauleap(algorithm.cache...)

  if isfinite(τ)
    rates = algorithm.rates
    # [OPTION 1] Use matrix-vector multiplication
    # 
    # Time: O(num_reactions + num_reactions * num_species)
    # Space: O(2*num_reactions + num_species)
    # 
    # map!(rate -> pois_rand(rate * τ), v, rates)
    # mul!(v, S, number_jumps)

    # [OPTION 2] Use explicit for-loop
    # 
    # Time: (mul/add) O(num_reactions * num_species) + (bool) O(num_reactions * num_species)
    # Space: O(num_species + num_reactions)
    # 
    # for j in eachindex(rates)
    #   num_jumps_j = pois_rand(rates[j] * τ)
    #   for k in eachindex(v)
    #     if j == 1
    #       v[k] = num_jumps_j * S[k, j]
    #     else
    #       v[k] += num_jumps_j * S[k, j]
    #     end
    #   end
    # end
  else
    fill!(v, zero(eltype(v)))
  end

  return v, s
end

function update!(algorithm::PoissonLeapAlgorithm, state, model, v)
  
end