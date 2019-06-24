# dirty hacks
get_number_jumps(model::ReactionSystem) = length(model.reactions)
get_number_jumps(model::InteractingParticleSystem) = length(model.reactions)
get_number_jumps(model::Vector{T}) where T = length(model)
get_number_species(state::Lattice) = length(state.types)
get_number_species(state::Vector{T}) where T = length(state)

abstract type SimulationAlgorithm end

struct Direct <: SimulationAlgorithm end

function build_simulator(::Direct, state, model, rates_cache)
  number_jumps = get_number_jumps(model)
  algorithm = DirectMethod{rates_cache}(zeros(number_jumps), 0.0)

  return ExactSimulator(algorithm)
end

struct EnhancedDirect <: SimulationAlgorithm end

function build_simulator(::EnhancedDirect, state, model, rates_cache)
  number_jumps = get_number_jumps(model)
  algorithm = EnhancedDirectMethod{HasRates}(zeros(number_jumps), 0.0)

  return ExactSimulator(algorithm)
end

struct SortingDirect <: SimulationAlgorithm end

function build_simulator(::SortingDirect, state, model, rates_cache)
  number_jumps = get_number_jumps(model)
  algorithm = SortingDirectMethod{HasRates}(zeros(number_jumps), 0.0, collect(1:number_jumps), 1)

  return ExactSimulator(algorithm)
end

struct FirstReaction <: SimulationAlgorithm end

function build_simulator(::FirstReaction, state, model, rates_cache)
  number_jumps = get_number_jumps(model)
  algorithm = FirstReactionMethod(zeros(number_jumps), 0.0)

  return ExactSimulator(algorithm)
end

struct NextReaction <: SimulationAlgorithm end

function build_simulator(::NextReaction, state, model, rates_cache)
  number_jumps = get_number_jumps(model)
  priority_queue = PQBinaryHeap{Int,Float64,ForwardOrdering}(collect(1:number_jumps), zeros(number_jumps))
  algorithm = NextReactionMethod(zeros(number_jumps), 0.0, priority_queue)

  return ExactSimulator(algorithm)
end

struct RejectionSSA <: SimulationAlgorithm end

function build_simulator(::RejectionSSA, state, model, rates_cache)
  number_species = length(state)
  number_jumps = get_number_jumps(model)

  interval_lo = Vector{Int}(undef, number_species)
  interval_hi = Vector{Int}(undef, number_species)
  rates = [Vector{Float64}(undef, 2) for j in 1:number_jumps]
  total_rate = Vector{Float64}(undef, 2)
  spc_dg = model.spc_graph
  rxn_dg = model.rxn_graph

  algorithm = RejectionMethod(state, model, interval_lo, interval_hi, rates, total_rate, spc_dg, rxn_dg)

  return ExactSimulator(algorithm)
end
