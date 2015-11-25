type SSA <: Algorithm
  itr::Int

  tf::Float64
  dt::Float64

  t::Float64
  intensity::Float64
  steps::Int

  function SSA(itr, tf, dt, args)
    new(itr, tf, dt, 0.0, 0.0, 0)
  end
end

init(alg::SSA, rrxns, spcs, initial, params) = return;

function reset(alg::SSA, rxns, spcs, params)
  alg.t     = 0.0
  alg.steps = 0
  return;
end

function step(alg::SSA, rxns, spcs, params)
  alg.intensity = compute_propensities!(rxns, spcs, params)
  τ = ssa_update!(spcs, rxns, alg.t, alg.tf, alg.intensity)
  alg.t = alg.t + τ
  alg.steps = alg.steps + 1
  return;
end

function sample(rxns::Vector{Reaction}, jump)
  ss = 0.0
  for i in eachindex(rxns)
    ss = ss + rxns[i].propensity
    if ss >= jump
      return i
    end
  end
  return 0
end

function ssa_step!(spcs::Vector{Int}, rxns::Vector{Reaction}, intensity)
  u = rand()
  jump = intensity * u
  j = sample(rxns, jump)
  j > 0 ? update!(spcs, rxns[j]) : error("No reaction occurred!")
  return;
end

function ssa_update!(spcs::Vector{Int}, rxns::Vector{Reaction}, t, tf, intensity)
  τ = rand(Exponential(1/intensity))
  t = t + τ
  if t > tf; return τ; end
  ssa_step!(spcs, rxns, intensity)
  return τ
end

function update!(spcs::Vector{Int}, r::Reaction)
  for i in eachindex(spcs)
    spcs[i] = spcs[i] + (r.post[i] - r.pre[i])
  end
  return;
end
