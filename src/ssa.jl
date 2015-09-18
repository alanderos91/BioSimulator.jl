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

function ssa(model::Simulation, t_final::Float64; itr::Int=1, tracing::Bool=false, output::Symbol=:explicit, stepsize::Float64=0.0)

  #Unpack model
  init = model.initial
  spcs = model.state
  rxns = model.rxns
  params = model.param

  md = Dict()

  # Results array
  results = Dict{ASCIIString, PopulationTrace}[]
  for i = 1:itr
    copy!(spcs, init)
    traces = init_traces(spcs)

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
    push!(results, traces)
  end

  md["algorithm"] = "SSA"
  md["iterations"] = itr
  md["duration"] = t_final
  if output == :explicit
    md["output"] = "Explicit"
  elseif output == :fixed && stepsize > 0.0
    md["output"] = "Fixed-Interval"
    md["stepsize"] = stepsize
    results = regularize(results, stepsize, t_final)
  end

  return SimulationResults(model.id, results, md)
end

function update!(spcs::Vector{Species}, r::Reaction)
  for i in 1:length(spcs)
    spcs[i].pop = spcs[i].pop + (r.post[i] - r.pre[i])
  end
  return;
end
