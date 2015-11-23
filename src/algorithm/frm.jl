function frm(model::Simulation, t_final::Float64, output::OutputType, dt::Float64, itr::Int)
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

    frm_steps = 0

    t = 0.0
    u.t_next = 0.0

    while t < t_final
      update!(output, df, t, u, spcs)
      compute_propensities!(rxns, spcs, params)
      τ = frm_update!(spcs, rxns, t, t_final)
      t = t + τ
      frm_steps = frm_steps + 1
    end
    final_update!(output, df, t, u, spcs)
  end
  return df
end

function frm_update!(spcs::Vector{Int}, rxns::Vector{Reaction}, t, t_final)
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
