export sal, dsal

function sal_step!(spcs::Vector{Species}, rxns::Vector{Reaction}, events::Vector{Int}, τ::Float64, α::Float64)
  while isbadleap(spcs, rxns, events)
    τ = contract!(events, τ, α)
  end

  update!(spcs, rxns, events)
  return τ
end

function sal(model::Simulation, t_final::Float64;
             itr::Int=1,
             tol::Float64=0.125,
             thrsh::Float64=100.0,
             ctrct::Float64=0.75)

  #Unpack model
  init = model.initial
  spcs = model.state
  rxns = model.rxns
  params = model.param

  job = SimulationJob(itr)

  # Create additional arrays
  drdt = zeros(Float64, length(rxns))
  events = zeros(Int, length(rxns))

  for i = 1:itr
    reset!(spcs, init)

    result = SimulationResult(spcs)
    ssa_steps = 0
    sal_steps = 0

    t = 0.0

    while t < t_final
      a_total = 0.0
      for r in rxns
        r.propensity = mass_action(r, spcs, params)
        a_total = a_total + r.propensity
      end

      τ = 0.0
      if a_total < thrsh
        τ = rand(Exponential(1/a_total))
        if t + τ <= t_final
          ssa_step!(spcs, rxns, a_total)
          ssa_steps = ssa_steps + 1
        end
      else
        compute_time_derivatives!(drdt, spcs, rxns, params)
        τ = tau_leap(rxns, params, drdt, tol)

        if t + τ <= t_final
          generate_events!(events, rxns, drdt, τ)
          τ = sal_step!(spcs, rxns, events, τ, ctrct)
          sal_steps = sal_steps + 1
        end
      end
      t = t + τ

      update!(result, t, spcs)
    end
    job[i] = result
  end
  return job
end

function dsal(model::Simulation, t_final::Float64;
             itr::Int=1,
             dt::Float64=1.0,
             tol::Float64=0.125,
             thrsh::Float64=100.0,
             ctrct::Float64=0.75)

  #Unpack model
  init = model.initial
  spcs = model.state
  rxns = model.rxns
  params = model.param

  n = round(Int, t_final / dt) + 1

  job = SimulationJob(itr)

  # Create additional arrays
  drdt = zeros(Float64, length(rxns))
  events = zeros(Int, length(rxns))

  for i = 1:itr
    reset!(spcs, init)

    result = SimulationResult(spcs, n)
    ssa_steps = 0
    sal_steps = 0

    t = 0.0
    t_next = dt
    j = 2

    while t < t_final
      a_total = 0.0
      for r in rxns
        r.propensity = mass_action(r, spcs, params)
        a_total = a_total + r.propensity
      end

      τ = 0.0
      if a_total < thrsh
        τ = rand(Exponential(1/a_total))
        if t + τ <= t_final
          ssa_step!(spcs, rxns, a_total)
          ssa_steps = ssa_steps + 1
        end
      else
        compute_time_derivatives!(drdt, spcs, rxns, params)
        τ = tau_leap(rxns, params, drdt, tol)

        if t + τ <= t_final
          generate_events!(events, rxns, drdt, τ)
          τ = sal_step!(spcs, rxns, events, τ, ctrct)
          sal_steps = sal_steps + 1
        end
      end
      t = t + τ

      while t >= t_next
        update!(result, t, spcs, j)
        j = j + 1
        t_next = t_next + dt
        if j > n; break; end
      end
    end
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

function compute_time_derivatives!(drdt::Vector{Float64}, spcs::Vector{Species}, rxns::Vector{Reaction}, param::Dict{ASCIIString, Float64})
  for i in eachindex(drdt)
    r = rxns[i]
    acc = 0.0
    for k in eachindex(spcs)
      ∂r∂x_k = mass_action_deriv(r, spcs, param, k)
      dxdt_k = mean_derivative(rxns, k)
      acc = acc + ∂r∂x_k * dxdt_k
    end
    drdt[i] = acc
  end
  return drdt
end

function mean_derivative(rxns::Vector{Reaction}, k::Int)
  acc = 0.0
  for j in eachindex(rxns)
    acc = acc + rxns[j].propensity * (rxns[j].post[k] - rxns[j].pre[k])
  end
  return acc
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
