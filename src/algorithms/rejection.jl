struct RejectionMethod{X,R,U,DG<:DependencyGraph} <: DecoupledSearchAlgorithm
  state::X
  model::R
  interval_lo::Vector{Int}
  interval_hi::Vector{Int}
  rates::Vector{U}
  total_rate::U
  spc_dg::DG
  rxn_dg::DG
end

# need to think about whether its possible to do this with cumulative rates
RatesCache(::RejectionMethod) = HasRates()

@inline function initialize!(algorithm::RejectionMethod, state, model)
  interval_lo = algorithm.interval_lo
  interval_hi = algorithm.interval_hi
  rates = algorithm.rates
  total_rate = algorithm.total_rate

  fill!(total_rate, zero(eltype(total_rate)))

  # build the upper and lower bounds on state
  for i in eachindex(state)
    interval_lo[i], interval_hi[i] = make_interval(state[i], 0.1)
  end

  # fill in the propensities

  for j in eachindex(rates)
    lo = rate(model, interval_lo, j)
    hi = rate(model, interval_hi, j)

    @inbounds rates[j][1] = lo
    @inbounds rates[j][2] = hi

    @inbounds total_rate[1] += lo
    @inbounds total_rate[2] += hi
  end
end

@inline function generate_jump(algorithm::RejectionMethod)
  j, logu = accept_reject(RatesCache(algorithm), algorithm.rates, algorithm.total_rate, algorithm.state, algorithm.model)
  # j, u = accept_reject(RatesCache(algorithm), algorithm.rates, algorithm.total_rate, algorithm.state, algorithm.model)
  s = logu / upper(algorithm.total_rate)
  # s = - log(u) / upper(algorithm.total_rate)

  return (j, s)
end

# rates should be an iterator over intervals_rates
@inline function accept_reject(ctype::S, rates, total_rate, state, model) where S
  logu = zero(eltype(total_rate))
  # r = one(eltype(total_rate))
  rejected = true
  j = 0

  while rejected
    # search jumps based on upper bounds
    @assert upper(total_rate) > 0 
    j = search_jump_rates(ctype, iterate_upper(rates), upper(total_rate))

    lo = lower(rates, j)
    hi = upper(rates, j)
    u  = rand() * hi

    rejected = (lo > zero(lo)) && (lo < u)

    if rejected
      rate_j = rate(model, state, j)
      rejected = (rate_j > zero(rate_j)) && (rate_j < u)
    end

    logu += randexp()
  end

  return j, logu
end

@inline function update!(algorithm::RejectionMethod, state, model, k)
  spc_dg = algorithm.spc_dg
  rxn_dg = algorithm.rxn_dg
  interval_lo = algorithm.interval_lo
  interval_hi = algorithm.interval_hi
  rates = algorithm.rates
  total_rate = algorithm.total_rate

  # for i in dependents(rxn_dg, k)
  for i in eachindex(state)
    @inbounds x = state[i]
    @inbounds x_lo = interval_lo[i]
    @inbounds x_hi = interval_hi[i]

    # check for states that have left interval
    if (x == 0) || (x < x_lo) || (x > x_hi)
      @inbounds interval_lo[i], interval_hi[i] = make_interval(x, 0.1)

      # check which propensity bounds must be updated based on a species-reaction dependency graph
      for j in dependents(spc_dg, i)
        @inbounds total_rate[1] -= rates[j][1]
        @inbounds total_rate[2] -= rates[j][2]

        @inbounds rates[j][1] = rate(model, interval_lo, j)
        @inbounds rates[j][2] = rate(model, interval_hi, j)

        @inbounds total_rate[1] += rates[j][1]
        @inbounds total_rate[2] += rates[j][2]
      end
    end
  end
end

function make_interval(x, p)
  if x == zero(x)
    lo = zero(x)
    hi = zero(x)
  elseif x < 25
    lo = max(0, x-4)
    hi = x+4
  else
    lo = trunc(Int, (1-p)*x)
    hi = trunc(Int, (1+p)*x)
  end

  return (lo, hi)
end

lower(x) = first(x)
upper(x) = last(x)

@inbounds lower(x, j) = first(x[j])
@inbounds upper(x, j) = last(x[j])

@inbounds iterate_upper(x) = (upper(x, j) for j in eachindex(x))
@inbounds iterate_lower(x) = (lower(x, j) for j in eachindex(x))
