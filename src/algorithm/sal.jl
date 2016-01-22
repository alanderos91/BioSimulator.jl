type SAL <: Algorithm
  tol::Float64
  thrsh::Float64
  ctrct::Float64

  intensity::Float64
  ssa_steps::Int
  sal_steps::Int

  dxdt::Vector{Float64}
  drdt::Vector{Float64}
  events::Vector{Int}

  function SAL(args)
    tol   = get(args, :tol,   0.125)
    thrsh = get(args, :thrsh, 100.0)
    ctrct = get(args, :ctrct, 0.75)

    new(tol, thrsh, ctrct, 0.0, 0, 0, Float64[], Float64[], Int[])
  end
end

function init(alg::SAL, rxns, spcs, initial, params)
  alg.dxdt   = zeros(Float64, length(spcs))
  alg.drdt   = zeros(Float64, length(rxns))
  alg.events = zeros(Int,     length(rxns))

  return;
end
function reset(alg::SAL, rxns, spcs, params)
  alg.ssa_steps = 0
  alg.sal_steps = 0
  return;
end

function step(alg::SAL, rxns, spcs, params, t, tf)
  alg.intensity = compute_propensities!(rxns, spcs, params)
  if alg.intensity < alg.thrsh
    τ = ssa_update!(spcs, rxns, t, tf, alg.intensity)
    alg.ssa_steps = alg.ssa_steps + 1
  else
    τ = sal_update!(spcs, rxns, t, tf, params, alg.dxdt, alg.drdt, alg.events, alg.tol, alg.ctrct)
    alg.sal_steps = alg.sal_steps + 1
  end
  return τ;
end

function update!(spcs::Vector{Int}, r::ReactionChannel, k::Int)
  for i in eachindex(spcs)
    spcs[i] = spcs[i] + k * (r.post[i] - r.pre[i])
  end
  return;
end

function update!(spcs::Vector{Int}, rxns::ReactionVector, events::Vector{Int})
  for i in eachindex(rxns)
    update!(spcs, rxns[i], events[i])
  end
  return;
end

function isbadleap(spcs::Vector{Int}, rxns::ReactionVector, events::Vector{Int})
  for i in eachindex(spcs)
    val = spcs[i]
    for j in eachindex(rxns)
      val = val + events[j] * (rxns[j].post[i] - rxns[j].pre[i])

      if val < 0
        return true
      end
    end
  end

  return false
end

function compute_time_derivatives!(drdt::Vector{Float64}, spcs::Vector{Int}, rxns::ReactionVector, param::Parameters, dxdt::Vector{Float64})
  for i in eachindex(drdt)
    drdt[i] = 0.0
    for k in eachindex(spcs)
      ∂r∂x_k = mass_action_deriv(rxns[i], spcs, param, k)
      drdt[i] = drdt[i] + ∂r∂x_k * dxdt[k]
    end
  end
  return drdt
end

function compute_mean_derivatives!(dxdt::Vector{Float64}, rxns::ReactionVector)
  for k in eachindex(dxdt)
    dxdt[k] = 0.0
    for j in eachindex(rxns)
      dxdt[k] = dxdt[k] + rxns[j].propensity * (rxns[j].post[k] - rxns[j].pre[k])
    end
  end
  return dxdt
end

function tau_leap(rxns::ReactionVector, param::Parameters, drdt::Vector{Float64}, ϵ::Float64)
  τ = Inf

  for j in eachindex(rxns)
    r = rxns[j]
    key = r.rate
    r_j = r.propensity
    c = param[key]

    a = ϵ * max(r_j, c)
    b = abs(drdt[j])

    if b > 0
      τ = min(τ, a / b)
    end
  end

  return τ
end

function generate_events!(events::Vector{Int}, rxns::ReactionVector, drdt::Vector{Float64}, τ::Float64)
  for i in eachindex(rxns)
    λ = τ * rxns[i].propensity + 0.5 * τ * τ * drdt[i]
    events[i] = rand(Poisson(λ))
  end
  return;
end

function contract!(events::Vector{Int}, τ::Float64, α::Float64)
  @inbounds for i in eachindex(events)
    k = 0
    @inbounds for j in 1:events[i]
      k = rand() < α ? k + 1 : k
    end
    events[i] = k
  end
  return τ * α
end

function sal_step!(spcs::Vector{Int}, rxns::ReactionVector, events::Vector{Int}, τ::Float64, α::Float64)
  while isbadleap(spcs, rxns, events)
    τ = contract!(events, τ, α)
  end

  update!(spcs, rxns, events)
  return τ
end

function sal_update!(spcs::Vector{Int}, rxns::ReactionVector, t, tf, params, dxdt, drdt, events, tol, ctrct)
  compute_mean_derivatives!(dxdt, rxns)
  compute_time_derivatives!(drdt, spcs, rxns, params, dxdt)
  τ = tau_leap(rxns, params, drdt, tol)
  τ = min(τ, tf - t)
  generate_events!(events, rxns, drdt, τ)
  τ = sal_step!(spcs, rxns, events, τ, ctrct)

  return τ
end
