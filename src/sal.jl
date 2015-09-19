export sal

function sal_step!(spcs::Vector{Species}, rxns::Vector{Reaction}, events::Vector{Int}, τ::Float64, α::Float64)
  while isbadleap(spcs, rxns, events)
    τ = contract!(events, τ, α)
  end

  update!(spcs, rxns, events)
  return τ
end

function sal(model::Simulation, t_final::Float64;
             itr::Int=1,
             tracing::Bool=false,
             output::Symbol=:explicit,
             stepsize::Float64=0.0,
             tol::Float64=0.125,
             thrsh::Float64=100.0,
             ctrct::Float64=0.75)

  #Unpack model
  init = model.initial
  spcs = model.state
  rxns = model.rxns
  params = model.param

  md = Dict()

  # Results array
  results = Dict{ASCIIString, PopulationTrace}[]

  # Create additional arrays
  drdt = zeros(Float64, length(rxns))
  events = zeros(Int, length(rxns))

  for i = 1:itr
    reset!(spcs, init)
    traces = init_traces(spcs)

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
    push!(results, traces)
  end

  set_metadata!(md, "SAL", t_final, itr, output, stepsize)
  if output == :fixed results = regularize(results, stepsize, t_final) end
  return SimulationResults(model.id, results, md)
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
  fill!(drdt, 0.0)

  for i in eachindex(drdt)
    for j in eachindex(rxns)
      r = rxns[j]
      c = param[rxns[j].rate]
      r_j = rxns[j].propensity

      for k in eachindex(r.pre)
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

  for j in eachindex(rxns)
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
  for i in eachindex(rxns)
    λ = τ * rxns[i].propensity + 0.5 * τ * τ * drdt[i]
    x = Poisson(λ)

    events[i] = round(Int, rand(x))
  end
  return;
end
