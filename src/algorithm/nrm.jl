function nrm(model::Simulation, t_final::Float64, output::OutputType, dt::Float64, itr::Int)
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

  g = init_dep_graph(rxns)
  pq = init_pq(rxns) # Allocate the priority queue with zeros

  for i = 1:itr
    copy!(spcs, init)

    nrm_steps = 0

    t = 0.0
    u.t_next = 0.0

    # Compute propensities and initialize the priority queue with firing times
    compute_propensities!(rxns, spcs, params)
    init_pq!(pq, rxns)

    while t < t_final
      update!(output, df, t, u, spcs)
      t = nrm_update!(spcs, rxns, t, t_final, g, pq, params)
      nrm_steps = nrm_steps + 1
    end
    final_update!(output, df, t, u, spcs)
  end
  return df
end

function nrm_update!(spcs::Vector{Int}, rxns::Vector{Reaction}, t::Float64, t_final::Float64, g::LightGraphs.DiGraph, pq::Collections.PriorityQueue{Int,Float64,Base.Order.ForwardOrdering}, param::Dict{ASCIIString,Float64})
  μ, t = Collections.peek(pq)
  if t > t_final; return t; end
  μ > 0 ? update!(spcs, rxns[μ]) : error("No reaction occurred!")
  update_dep_graph!(g, rxns, pq, spcs, param, μ, t)
  return t
end

function init_dep_graph(rxns::Vector{Reaction})
  m = length(rxns)
  d = length(rxns[1].pre)
  g = DiGraph(m)

  @inbounds for j in eachindex(rxns)
    r = rxns[j]
    pre = r.pre
    post = r.post

    # Search for reactions r_i with reactant or product of r_j as reactant
    @inbounds for k in eachindex(pre)
      if pre[k] != 0 || post[k] != 0
        @inbounds for i in eachindex(rxns)
          if rxns[i].pre[k] != 0; add_edge!(g, j, i); end
        end
      end
    end
  end
  return g
end

function update_dep_graph!(g::LightGraphs.DiGraph, rxns::Vector{Reaction}, pq::Collections.PriorityQueue{Int,Float64,Base.Order.ForwardOrdering}, spcs::Vector{Int}, param::Dict{ASCIIString,Float64}, μ::Int, t::Float64)
  dependents = neighbors(g, μ)
  r = rxns[μ]

  # Compute propensity of the reaction that fired and update the next firing time
  propensity!(r, spcs, param)
  pq[μ] = t + rand(Exponential(1 / r.propensity))

  # Compute propensities and firing times for dependent reactions
  @inbounds for α in dependents
    _r = rxns[α]
    old_propensity = _r.propensity
    propensity!(_r, spcs, param)
    if α != μ
      pq[α] = t + (old_propensity / _r.propensity) * (pq[α] - t)
    end
  end
  return g
end

function init_pq(rxns::Vector{Reaction})
  pq = Collections.PriorityQueue(Int,Float64)
  for j in eachindex(rxns)
    pq[j] = 0.0
  end
  return pq
end

function init_pq!(pq::Collections.PriorityQueue{Int,Float64,Base.Order.ForwardOrdering}, rxns::Vector{Reaction})
  @inbounds for j in eachindex(rxns)
    pq[j] = rand(Exponential(1 / rxns[j].propensity))
  end
end
