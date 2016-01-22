type FRM <: Algorithm
  intensity::Float64
  steps::Int

  function FRM(args)
    new(0.0, 0)
  end
end

init(alg::FRM, rxns, spcs, initial, params) = return;

function reset(alg::FRM, rxns, spcs, params)
  alg.steps = 0

  return;
end

function step(alg::FRM, rxns, spcs, params, t, tf)
  compute_propensities!(rxns, spcs, params)
  τ = frm_update!(spcs, rxns, t, tf)
  alg.steps = alg.steps + 1

  return τ;
end

function frm_update!(spcs::Vector{Int}, rxns::ReactionVector, t, tf)
  τ = Inf; μ = 0
  for j in eachindex(rxns)
    τj = rand(Exponential(1/rxns[j].propensity))
    if τj < τ
      τ = τj
      μ = j
    end
  end
  t = t + τ
  if t > tf; return τ; end
  μ > 0 ? update!(spcs, rxns[μ]) : error("No reaction occurred!")
  return τ
end
