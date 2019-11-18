mutable struct BinomialLeapMethod{F1,U} <: SafeLeapAlgorithm
  rates::Vector{Float64}
  total_rate::Float64
  stoichiometry::SparseMatrixCSC{Int,Int}
  leap_formula::F1
  available_reactants::U
end

# information we need:
# - binomial distributions
#     - maximum number of firings for reaction j
#     - sampling probability for reaction j
# - available reactants
# - leap size
#     - tau = max_j (max # of firings for reaction j) / rate_j
