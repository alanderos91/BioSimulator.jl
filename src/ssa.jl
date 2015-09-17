function ssa_step!(spcs::Vector{Species}, rxns::Vector{Reaction}, jump::Float64)
  ss = 0.0
  flag = true

  for r in rxns
    ss = ss + r.propensity
    if jump <= ss
      update!(spcs, r)
      flag = false
      break
    end
  end

  if flag
    error("No reaction ocurred...")
  end
end

function ssa!(spcs::Vector{Species}, rxns::Vector{Reaction}, params::Dict{ASCIIString, Float64}, t_final::Float64)
  t = 0.0

  while t < t_final
    u1 = rand()
    u2 = rand()
    a_total = 0.0

    for r in rxns
      r.propensity = mass_action(r, spcs, params)
      a_total = a_total + r.propensity
    end

    τ = -log(u1) / a_total
    jump = a_total * u2

    ssa_step!(spcs, rxns, jump)
    t = t + τ
  end

  return spcs
end

function ssa_explicit(model::Simulation, t_final::Float64; tracing::Bool=false)

  traces = Dict{ASCIIString, PopulationTrace}()
  for s in model.initial
    if s.istracked
      traces[s.id] = PopulationTrace(s)
    end
  end

  #Unpack model
  spcs = model.state
  rxns = model.rxns
  params = model.param

  t = 0.0

  while t < t_final
    u1 = rand()
    u2 = rand()
    a_total = 0.0

    for r in rxns
      r.propensity = mass_action(r, spcs, params)
      a_total = a_total + r.propensity
    end

    τ = -log(u1) / a_total
    if t + τ <= t_final
      jump = a_total * u2
      ssa_step!(spcs, rxns, jump)
    end

    t = t + τ

    update_traces!(traces, t, spcs, tracing)
  end

  result = SimulationResult("SSA (Explicit)", traces, Dict{Any,Any}())

  result.metadata["duration"] = t_final

  return result
end

function ssa_stepper!(model::Simulation, steps::Int, stepsize::Float64; tracing::Bool=false)
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

  for i = 1:steps
    ssa!(spcs, rxns, params, stepsize)

    update_traces!(traces, i * stepsize, spcs, tracing)
  end

  result = SimulationResult("SSA (Fixed)", traces, Dict{Any,Any}())

  result.metadata["duration"] = steps * stepsize
  result.metadata["stepsize"] = stepsize

  return result
end

function update!(spcs::Vector{Species}, r::Reaction)
  for i in 1:length(spcs)
    spcs[i].pop = spcs[i].pop + (r.post[i] - r.pre[i])
  end
  return;
end
