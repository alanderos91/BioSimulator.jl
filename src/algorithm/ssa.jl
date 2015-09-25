export ssa, dssa

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

function ssa(model::Simulation, t_final::Float64; itr::Int=1)

  #Unpack model
  init = model.initial
  spcs = model.state
  rxns = model.rxns
  params = model.param

  job = SimulationJob(itr)

  for i = 1:itr
    reset!(spcs, init)

    result = SimulationResult(spcs)
    ssa_steps = 0

    t = 0.0

    while t < t_final
      a_total = 0.0
      for r in rxns
        propensity!(r, spcs, params)
        a_total = a_total + r.propensity
      end

      τ = rand(Exponential(1/a_total))
      if t + τ <= t_final
        ssa_step!(spcs, rxns, a_total)
        ssa_steps = ssa_steps + 1
      end
      t = t + τ

      update!(result, t, spcs)
    end
    job[i] = result
  end
  return job
end

function dssa(model::Simulation, t_final::Float64; itr::Int=1, dt::Float64=1.0)
  #Unpack model
  init = model.initial
  spcs = model.state
  rxns = model.rxns
  params = model.param

  n = round(Int, t_final / dt) + 1

  job = SimulationJob(itr)

  for i = 1:itr
    reset!(spcs, init)

    result = SimulationResult(spcs, n)
    ssa_steps = 0

    t = 0.0
    t_next = dt
    j = 2

    while t < t_final
      a_total = 0.0
      for r in rxns
        propensity!(r, spcs, params)
        a_total = a_total + r.propensity
      end

      τ = rand(Exponential(1/a_total))
      if t + τ <= t_final
        ssa_step!(spcs, rxns, a_total)
        ssa_steps = ssa_steps + 1
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

function update!(spcs::Vector{Species}, r::Reaction)
  for i in eachindex(spcs)
    spcs[i].pop = spcs[i].pop + (r.post[i] - r.pre[i])
  end
  return;
end
