export nrm, dnrm

function nrm(model::Simulation, t_final::Float64; itr::Int=1)
  #Unpack model
  init = model.initial
  spcs = model.state
  rxns = model.rxns
  params = model.param

  job = SimulationJob(itr)

  g = init_dep_graph(rxns)
  pq = init_pq(rxns) # Allocate the priority queue with zeros

  for i = 1:itr
    reset!(spcs, init)

    result = SimulationResult(spcs)
    nrm_steps = 0

    t = 0.0

    # Compute propensities and initialize the priority queue with firing times
    compute_propensities!(rxns, spcs, params)
    init_pq!(pq, rxns)

    while t < t_final
      update!(result, t, spcs)
      t = nrm_update!(spcs, rxns, t, t_final, g, pq, params)
      nrm_steps = nrm_steps + 1
    end
    update!(result, t_final, spcs)
    job[i] = result
  end
  return job
end

function dnrm(model::Simulation, t_final::Float64; itr::Int=1, dt::Float64=1.0)
  #Unpack model
  init = model.initial
  spcs = model.state
  rxns = model.rxns
  params = model.param

  n = round(Int, t_final / dt) + 1

  job = SimulationJob(itr)

  g = init_dep_graph(rxns)
  pq = init_pq(rxns) # Allocate the priority queue with zeros

  for i = 1:itr
    reset!(spcs, init)

    result = SimulationResult(spcs, n)
    nrm_steps = 0

    t = 0.0
    t_next = 0.0
    j = 1

    # Compute propensities and initialize the priority queue with firing times
    compute_propensities!(rxns, spcs, params)
    init_pq!(pq, rxns)

    while t < t_final
      while t >= t_next
        update!(result, t_next, spcs, j)
        j = j + 1
        t_next = t_next + dt
        if j > n; break; end
      end
      t = nrm_update!(spcs, rxns, t, t_final, g, pq, params)
      nrm_steps = nrm_steps + 1
    end
    while j <= n
      update!(result, t_next, spcs, j)
      j = j + 1
      t_next = t_next + dt
    end
    job[i] = result
  end
  return job
end

function nrm_update!(spcs::Vector{Species}, rxns::Vector{Reaction}, t, t_final, g, pq, param)
  μ, t = Collections.peek(pq)
  if t > t_final; return t; end
  μ > 0 ? update!(spcs, rxns[μ]) : error("No reaction occurred!")
  update_dep_graph!(g, rxns, pq, spcs, param, μ, t)
  return t
end

function init_dep_graph(rxns::Vector{Reaction})
  g = DiGraph(length(rxns))
  for j in eachindex(rxns)
    r = rxns[j]
    #add_edge!(g, j, j)

    # Search for reactions R_i with reactant or product of R_j as reactant
    for k in eachindex(r.pre)
      if (r.pre[k] == 0) && (r.post[k] == 0); continue; end

      for i in eachindex(rxns)
        #if i == j; continue; end
        if (rxns[i].pre[k] != 0); add_edge!(g, j, i); end
      end
    end
  end
  return g
end

function update_dep_graph!(g, rxns::Vector{Reaction}, pq, spcs::Vector{Species}, param, μ, t)
  dependents = neighbors(g, μ)

  # Compute propensity of the reaction that fired and update the next firing time
  propensity!(rxns[μ], spcs, param)
  pq[μ] = t + rand(Exponential(1 / rxns[μ].propensity))

  # Compute propensities and firing times for dependent reactions
  for α in dependents
    old_propensity = rxns[α].propensity
    propensity!(rxns[α], spcs, param)
    if α != μ
      pq[α] = t + (old_propensity / rxns[α].propensity) * (pq[α] - t)
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

function init_pq!(pq, rxns::Vector{Reaction})
  for j in eachindex(rxns)
    pq[j] = rand(Exponential(1 / rxns[j].propensity))
  end
end
