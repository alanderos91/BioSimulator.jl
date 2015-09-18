export ssa

function sample(rxns::Vector{Reaction}, jump::Float64)
  ss = 0.0
  for i = 1:length(rxns)
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

function ssa(model::Simulation, t_final::Float64; tracing::Bool=false, output::Symbol=:explicit, stepsize::Float64=0.0)

  traces = Dict{ASCIIString, PopulationTrace}()
  for s in model.initial
    if s.istracked
      traces[s.id] = PopulationTrace(s)
    end
  end
  md = Dict()

  #Unpack model
  spcs = model.state
  rxns = model.rxns
  params = model.param

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
    end
    t = t + τ

    update_traces!(traces, t, spcs, tracing)
  end

  md["duration"] = t_final
  if output == :explicit
    md["output"] = "Explicit"
  elseif output == :fixed && stepsize > 0.0
    md["output"] = "Fixed-Interval"
    md["stepsize"] = stepsize
    traces = regularize(traces, stepsize, t_final)
  end

  return SimulationResult("SSA", traces, md)
end

function update!(spcs::Vector{Species}, r::Reaction)
  for i in 1:length(spcs)
    spcs[i].pop = spcs[i].pop + (r.post[i] - r.pre[i])
  end
  return;
end
