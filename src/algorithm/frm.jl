export frm

function frm(model::Simulation, t_final::Float64; itr::Int=1)
  #Unpack model
  init = model.initial
  spcs = model.state
  rxns = model.rxns
  params = model.param

  job = SimulationJob(itr)

  for i = 1:itr
    reset!(spcs, init)

    result = SimulationResult(spcs)
    frm_steps = 0

    t = 0.0

    while t < t_final
      update!(result, t, spcs)
      compute_propensities!(rxns, spcs, params)
      τ = frm_update!(spcs, rxns, t, t_final)
      t = t + τ
      frm_steps = frm_steps + 1
    end
    update!(result, t_final, spcs)
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
