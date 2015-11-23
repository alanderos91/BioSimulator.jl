function ssa(model::Simulation, t_final::Float64, output::OutputType, dt::Float64, itr::Int)
  # Unpack model
  init = model.initial
  sname = model.sname
  tracked = model.tracked
  spcs = model.state
  rxns = model.rxns
  params = model.param

  n = round(Int, t_final / dt) + 1

  u  = if isa(output, Explicit)
         Updater(dt, tracked, 0)
       elseif isa(output, Uniform)
         Updater(dt, tracked, n)
       elseif isa(output, Mean)
         Updater(dt, tracked, n)
       elseif isa(output, Histogram)
         Updater(dt, tracked, itr)
       end

  df = if isa(output, Explicit)
         init_df(sname, tracked, 0)
       elseif isa(output, Uniform)
         init_df(sname, tracked, itr * n)
       elseif isa(output, Mean)
         init_df(sname, tracked, n)
       elseif isa(output, Histogram)
         init_df(sname, tracked, itr)
       end

  for i = 1:itr
    copy!(spcs, init)

    ssa_steps = 0

    t = 0.0

    while t < t_final
      update!(output, df, t, u, spcs)
      intensity = compute_propensities!(rxns, spcs, params)
      τ = ssa_update!(spcs, rxns, t, t_final, intensity)
      t = t + τ
      ssa_steps = ssa_steps + 1
    end
    final_update!(output, df, t, u, spcs)
  end
  return df
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

function ssa_step!(spcs::Vector{Int}, rxns::Vector{Reaction}, a_total::Float64)
  u = rand()
  jump = a_total * u
  j = sample(rxns, jump)
  j > 0 ? update!(spcs, rxns[j]) : error("No reaction occurred!")
  return;
end

function ssa_update!(spcs::Vector{Int}, rxns::Vector{Reaction}, t, t_final, intensity)
  τ = rand(Exponential(1/intensity))
  t = t + τ
  if t > t_final; return τ; end
  ssa_step!(spcs, rxns, intensity)
  return τ
end

function update!(spcs::Vector{Int}, r::Reaction)
  for i in eachindex(spcs)
    spcs[i] = spcs[i] + (r.post[i] - r.pre[i])
  end
  return;
end
