function sal_step!(spcs::Vector{Species}, rxns::Vector{Reaction}, events::Vector{Int}, τ::Float64, α::Float64)
  while isbadleap(spcs, rxns, events)
    τ = contract!(events, τ, α)
  end

  update!(spcs, rxns, events)
  return τ
end

function sal!(spcs::Vector{Species},
              rxns::Vector{Reaction},
              param::Dict{ASCIIString, Float64},
              t_final::Float64,
              ϵ::Float64,
              δ::Float64,
              α::Float64,
              drdt::Vector{Float64},
              events::Vector{Int})
  t = 0.0

  while t < t_final

    a_total = 0.0
    for r in rxns
      r.propensity = mass_action(r, spcs, param)
      a_total = a_total + r.propensity
    end

    τ = 0.0
    if a_total < δ
      u1 = rand()
      u2 = rand()
      τ = -log(u1) / a_total
      jump = a_total * u2
      ssa_step!(spcs, rxns, jump)
    else
      compute_time_derivatives!(drdt, spcs, rxns, param)
      τ = tau_leap(rxns, param, drdt, ϵ)

      generate_events!(events, rxns, drdt, τ)
      τ = sal_step!(spcs, rxns, events, τ, α)
    end

    t = t + τ
  end

  return spcs
end

function sal_explicit(model::Simulation, t_final::Float64;
                      tracing::Bool=false,
                      tol::Float64=0.125,
                      thrsh::Float64=100.0,
                      ctrct::Float64=0.75)
  traces = Dict{ASCIIString, PopulationTrace}()
  for s in model.initial
    if s.istracked
      traces[s.id] = PopulationTrace(s)
    end
  end

  # Unpack model
  spcs = model.state
  rxns = model.rxns
  params = model.param

  # Create additional arrays
  drdt = zeros(Float64, length(rxns))
  events = zeros(Int, length(rxns))

  t = 0.0
  while t < t_final
    a_total = 0.0
    for r in rxns
      r.propensity = mass_action(r, spcs, params)
      a_total = a_total + r.propensity
    end

    τ = 0.0
    if a_total < thrsh
      u1 = rand()
      u2 = rand()
      τ = -log(u1) / a_total
      if t + τ <= t_final
        jump = a_total * u2
        ssa_step!(spcs, rxns, jump)
      end
    else
      compute_time_derivatives!(drdt, spcs, rxns, params)
      τ = tau_leap(rxns, params, drdt, tol)

      if t + τ <= t_final
        generate_events!(events, rxns, drdt, τ)
        τ = sal_step!(spcs, rxns, events, τ, ctrct)
      end
    end

    t = t + τ
    update_traces!(traces, t, spcs, tracing)
  end

  result = SimulationResult("SAL (Fixed)", traces, Dict{Any,Any}())

  result.metadata["duration"] = t_final

  return result
end

function sal_stepper!(model::Simulation, steps::Int, stepsize::Float64;
                      tracing::Bool=false,
                      tol::Float64=0.125,
                      thrsh::Float64=100.0,
                      ctrct::Float64=0.75)
  # TODO: How to initialize wrt to `tracing` var

  traces = Dict{ASCIIString, PopulationTrace}()
  for s in model.initial
    if s.istracked
      traces[s.id] = PopulationTrace(s)
    end
  end

  # Unpack model
  spcs = model.state
  rxns = model.rxns
  params = model.param

  # Create additional arrays
  drdt = zeros(Float64, length(rxns))
  events = zeros(Int, length(rxns))

  for i = 1:steps
    sal!(spcs, rxns, params, stepsize, tol, thrsh, ctrct, drdt, events)

    update_traces!(traces, i * stepsize, spcs, tracing)
  end

  result = SimulationResult("SAL (Fixed)", traces, Dict{Any,Any}())

  result.metadata["duration"] = steps * stepsize
  result.metadata["stepsize"] = stepsize

  return result
end

function update!(spcs::Vector{Species}, r::Reaction, k::Int)
  for i = 1:length(spcs)
    spcs[i].pop = spcs[i].pop + k * (r.post[i] - r.pre[i])
  end
  return;
end

function update!(spcs::Vector{Species}, rxns::Vector{Reaction}, events::Vector{Int})
  for i = 1:length(rxns)
    update!(spcs, rxns[i], events[i])
  end
  return;
end

function isbadleap(spcs::Vector{Species}, rxns::Vector{Reaction}, events::Vector{Int})
  for i = 1:length(spcs)
    val = spcs[i].pop
    for j = 1:length(rxns)
      val = val + events[j] * (rxns[j].post[i] - rxns[j].pre[i])

      if val < 0
        return true
      end
    end
  end

  return false
end

function compute_time_derivatives!(drdt::Vector{Float64}, spcs::Vector{Species}, rxns::Vector{Reaction}, param::Dict{ASCIIString, Float64})
  fill!(drdt, 0.0)

  for i = 1:length(drdt)
    for j = 1:length(rxns)
      r = rxns[j]
      c = param[rxns[j].rate]
      r_j = rxns[j].propensity

      for k = 1:length(r.pre)
        ν = rxns[j].post[k] - rxns[j].pre[k]
        ∂r∂x = mass_action_deriv(c, r.pre, spcs, k)
        drdt[i] = drdt[i] + r_j * ∂r∂x * ν
      end
    end
  end
  return;
end

function tau_leap(rxns::Vector{Reaction}, param::Dict{ASCIIString, Float64}, drdt::Vector{Float64}, ϵ::Float64)
  τ = Inf
  n = length(rxns)

  for j = 1:n
    r = rxns[j]
    ratesym::ASCIIString = r.rate
    r_j = r.propensity
    c::Float64 = param[ratesym]

    a = ϵ * max(r_j, c)
    b = abs(drdt[j])

    if b > 0
      τ = min(τ, a / b)
    end
  end

  return τ
end

function generate_events!(events::Vector{Int}, rxns::Vector{Reaction}, drdt::Vector{Float64}, τ::Float64)
  for i = 1:length(rxns)
    λ = τ * rxns[i].propensity + 0.5 * τ * τ * drdt[i]
    x = Poisson(λ)

    events[i] = round(Int, rand(x))
  end
  return;
end
