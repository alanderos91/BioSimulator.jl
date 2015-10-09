function ssa(model::Simulation, t_final::Float64, output::OutputType, dt::Float64, itr::Int)
  # Unpack model
  init = model.initial
  spcs = model.state
  rxns = model.rxns
  params = model.param

  n = round(Int, t_final / dt) + 1

  job = SimulationJob(itr)

  for i = 1:itr
    reset!(spcs, init)

    result = init_sr(output, spcs, n)
    ssa_steps = 0

    t = 0.0
    t_next = 0.0
    j = 1

    while t < t_final
      t_next, j = update!(output, result, n, t, t_next, dt, j, spcs)
      intensity = compute_propensities!(rxns, spcs, params)
      τ = ssa_update!(spcs, rxns, t, t_final, intensity)
      t = t + τ
      ssa_steps = ssa_steps + 1
    end
    update!(output, result, n, t, t_final, dt, j, spcs)
    job[i] = result
  end
  return job
end

function sample(rxns::Vector{Reaction}, jump::Float64)
  ss = 0.0
  for i in eachindex(rxns)
    ss = ss + rxns[i].propensity
    if ss >= jump
      return i
    end
  end
  return 0
end

function ssa_step!(spcs::Vector{Species}, rxns::Vector{Reaction}, a_total::Float64)
  u = rand()
  jump = a_total * u
  j = sample(rxns, jump)
  j > 0 ? update!(spcs, rxns[j]) : error("No reaction occurred!")
  return;
end

function ssa_update!(spcs::Vector{Species}, rxns::Vector{Reaction}, t, t_final, intensity)
  τ = rand(Exponential(1/intensity))
  t = t + τ
  if t > t_final; return τ; end
  ssa_step!(spcs, rxns, intensity)
  return τ
end

function update!(spcs::Vector{Species}, r::Reaction)
  for i in eachindex(spcs)
    spcs[i].pop = spcs[i].pop + (r.post[i] - r.pre[i])
  end
  return;
end
