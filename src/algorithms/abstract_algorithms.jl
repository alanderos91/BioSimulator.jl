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

@inline function update_jump_rates!(algorithm::UnsafeLeapAlgorithm, state, model)
  algorithm.total_rate = update_jump_rates!(RatesCache(algorithm), algorithm.rates, algorithm.total_rate, state, model)
end

function initialize!(algorithm::UnsafeLeapAlgorithm, state, model)
  update!(algorithm, state, model)
end

function generate_leap!(v, algorithm::UnsafeLeapAlgorithm, Δt_max)
  # unpack data
  rates = algorithm.rates
  total_rate = algorithm.total_rate
  validate_leap! = algorithm.validate_leap!

  # evaluate the leap formula to compute τ
  s = algorithm.leap_formula(rates, total_rate)
  s = min(s, Δt_max)

  # if there is a valid leap...
  if isfinite(s)
    generate_events!(algorithm, v, rates, s)
    validate_leap!(v, s)
  else
    # otherwise fill the jumps with zeros
    fill!(v, zero(eltype(v)))
  end

  return v, s
end
