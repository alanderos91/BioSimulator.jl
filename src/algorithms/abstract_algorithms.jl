abstract type AbstractAlgorithm end

# trait dispatch
RatesCache(::AbstractAlgorithm) = HasRates()

# by default, initialization is the same as updating the rates vector
# algorithms should overload this behavior as needed
@inline function initialize!(algorithm::AbstractAlgorithm, state, model)
  update_jump_rates!(algorithm, state, model)
end

# we assume that `s` is a time increment rather than an absolute time
# as such, the default behavior is to return `t + s` as the next system time
# algorithms such as NRM can overload this behavior
@inline @fastmath function get_new_time(algorithm::AbstractAlgorithm, t, s)
  return t + s
end

abstract type DecoupledSearchAlgorithm <: AbstractAlgorithm end

# delegate the update procedure based on the RatesCacheType
@inline function update_jump_rates!(algorithm::DecoupledSearchAlgorithm, state, model)
  algorithm.total_rate = update_jump_rates!(RatesCache(algorithm), algorithm.rates, algorithm.total_rate, state, model)
end

@inline function update_jump_rates!(algorithm::DecoupledSearchAlgorithm, state, model, j)
  algorithm.total_rate = update_jump_rates!(RatesCache(algorithm), algorithm.rates, algorithm.total_rate, state, model, j)
end

# A DecoupledSearchSampler always generates the next reaction time based on an exponential distribution.
# The next reaction index depends exclusively on the sampler.
@inline function generate_jump(algorithm::DecoupledSearchAlgorithm)
  j = search_jumps(algorithm)
  s = randexp() / algorithm.total_rate

  return (j, s)
end

abstract type CoupledSearchAlgorithm <: AbstractAlgorithm end

# delegate the update procedure based on the RatesCacheType
@inline function update_jump_rates!(algorithm::CoupledSearchAlgorithm, state, model)
  algorithm.total_rate = update_jump_rates!(RatesCache(algorithm), algorithm.rates, algorithm.total_rate, state, model)
end

@inline function update_jump_rates!(algorithm::CoupledSearchAlgorithm, state, model, j)
  algorithm.total_rate = update_rates!(RatesCache(algorithm), algorithm.rates, algorithm.total_rate, state, model, j)
end

@inline function generate_jump(algorithm::CoupledSearchAlgorithm)
  (j, s) = search_jumps(algorithm)

  return (j, s)
end

abstract type UnsafeLeapAlgorithm <: AbstractAlgorithm end

##### recovering from negative excursions

function generate_leap!(v, algorithm::UnsafeLeapAlgorithm)
  # unpack data
  rates = algorithm.rates
  tota_rate = algorithm.total_rate

  # evaluate the leap formula to compute τ
  s = algorithm.leap_formula(rates, total_rate)

  # if there is a valid leap...
  if isfinite(τ)
    generate_events!(algorithm, v, rates, s)
  else
    # otherwise fill the jumps with zeros
    fill!(v, zero(eltype(v)))
  end

  return v, s
end

@inline function validate_leap!(simulator::TauLeapSimulator, state, model)
  Δ = simulator.proposed_change
  A = simulator.stoichiometry
  v = simulator.next_leap_jumps
  s = simulator.next_leap_time

  # check and contract invalid leaps
  mul!(Δ, A, v)
  invalid_leap = is_invalid_leap(Δ, state)

  while invalid_leap
    algorithm.contract!(v, s)
    mul!(proposal, stoichiometry, v)
    invalid_leap = is_invalid_leap(proposal, state)
  end
end