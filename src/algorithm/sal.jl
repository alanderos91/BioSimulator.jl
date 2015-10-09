function sal(model::Simulation, t_final::Float64, output::OutputType, dt::Float64, itr::Int, tol::Float64, thrsh::Float64, ctrct::Float64)
  #Unpack model
  init = model.initial
  spcs = model.state
  rxns = model.rxns
  params = model.param

  n = round(Int, t_final / dt) + 1

  job = SimulationJob(itr)

  # Create additional arrays
  dxdt = zeros(Float64, length(spcs)) # Rates of change for mean particle counts
  drdt = zeros(Float64, length(rxns)) # Rates of change for reaction propensities
  events = zeros(Int, length(rxns))   # Number of reaction events for each channel

  for i = 1:itr
    reset!(spcs, init)

    result = init_sr(output, spcs, n)
    ssa_steps = 0
    sal_steps = 0

    t = 0.0
    t_next = 0.0
    j = 1

    while t < t_final
      t_next, j = update!(output, result, n, t, t_next, dt, j, spcs)
      intensity = compute_propensities!(rxns, spcs, params)

      if intensity < thrsh
        τ = ssa_update!(spcs, rxns, t, t_final, intensity)
        t = t + τ
        ssa_steps = ssa_steps + 1
      else
        τ = sal_update!(spcs, rxns, t, t_final, params, dxdt, drdt, events, tol, ctrct)
        t = t + τ
        sal_steps = sal_steps + 1
      end
    end
    update!(output, result, n , t, t_final, dt, j, spcs)
    job[i] = result
  end
  return job
end

function update!(spcs::Vector{Species}, r::Reaction, k::Int)
  for i in eachindex(spcs)
    spcs[i].pop = spcs[i].pop + k * (r.post[i] - r.pre[i])
  end
  return;
end

function update!(spcs::Vector{Species}, rxns::Vector{Reaction}, events::Vector{Int})
  for i in eachindex(rxns)
    update!(spcs, rxns[i], events[i])
  end
  return;
end

function isbadleap(spcs::Vector{Species}, rxns::Vector{Reaction}, events::Vector{Int})
  for i in eachindex(spcs)
    val = spcs[i].pop
    for j in eachindex(rxns)
      val = val + events[j] * (rxns[j].post[i] - rxns[j].pre[i])

      if val < 0
        return true
      end
    end
  end

  return false
end

function compute_time_derivatives!(drdt::Vector{Float64}, spcs::Vector{Species}, rxns::Vector{Reaction}, param::Dict{ASCIIString, Float64}, dxdt::Vector{Float64})
  for i in eachindex(drdt)
    drdt[i] = 0.0
    for k in eachindex(spcs)
      ∂r∂x_k = mass_action_deriv(rxns[i], spcs, param, k)
      drdt[i] = drdt[i] + ∂r∂x_k * dxdt[k]
    end
  end
  return drdt
end

function compute_mean_derivatives!(dxdt::Vector{Float64}, rxns::Vector{Reaction})
  for k in eachindex(dxdt)
    dxdt[k] = 0.0
    for j in eachindex(rxns)
      dxdt[k] = dxdt[k] + rxns[j].propensity * (rxns[j].post[k] - rxns[j].pre[k])
    end
  end
  return dxdt
end

function tau_leap(rxns::Vector{Reaction}, param::Dict{ASCIIString, Float64}, drdt::Vector{Float64}, ϵ::Float64)
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

function generate_events!(events::Vector{Int}, rxns::Vector{Reaction}, drdt::Vector{Float64}, τ::Float64)
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

function sal_step!(spcs::Vector{Species}, rxns::Vector{Reaction}, events::Vector{Int}, τ::Float64, α::Float64)
  while isbadleap(spcs, rxns, events)
    τ = contract!(events, τ, α)
  end

  update!(spcs, rxns, events)
  return τ
end

function sal_update!(spcs::Vector{Species}, rxns::Vector{Reaction}, t, t_final, params, dxdt, drdt, events, tol, ctrct)
  compute_mean_derivatives!(dxdt, rxns)
  compute_time_derivatives!(drdt, spcs, rxns, params, dxdt)
  τ = tau_leap(rxns, params, drdt, tol)
  τ = min(τ, t_final - t)
  generate_events!(events, rxns, drdt, τ)
  τ = sal_step!(spcs, rxns, events, τ, ctrct)

  return τ
end
