function frm(model::Simulation, t_final::Float64, output::OutputType, dt::Float64, itr::Int)
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
    frm_steps = 0

    t = 0.0
    t_next = 0.0
    j = 1

    while t < t_final
      t_next, j = update!(output, result, n, t, t_next, dt, j, spcs)
      compute_propensities!(rxns, spcs, params)
      τ = frm_update!(spcs, rxns, t, t_final)
      t = t + τ
      frm_steps = frm_steps + 1
    end
    update!(output, result, n, t_next, dt, j, spcs)
    job[i] = result
  end
  return job
end

function frm_update!(spcs::Vector{Species}, rxns::Vector{Reaction}, t, t_final)
  τ = Inf; μ = 0
  for j in eachindex(rxns)
    τj = rand(Exponential(1/rxns[j].propensity))
    if τj < τ
      τ = τj
      μ = j
    end
  end
  t = t + τ
  if t > t_final; return τ; end
  μ > 0 ? update!(spcs, rxns[μ]) : error("No reaction occurred!")
  return τ
end
