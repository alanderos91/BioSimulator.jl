abstract type SimulationAlgorithm end

struct Direct <: SimulationAlgorithm end

function build_simulator(::Direct, state, model, rates_cache)
  number_jumps = length(model.reactions)
  algorithm = DirectMethod{rates_cache}(zeros(number_jumps), 0.0)

  return ExactSimulator(algorithm)
end

struct EnhancedDirect <: SimulationAlgorithm end

function build_simulator(::EnhancedDirect, state, model, rates_cache)
  number_jumps = length(model.reactions)
  algorithm = EnhancedDirectMethod{HasRates}(zeros(number_jumps), 0.0)

  return ExactSimulator(algorithm)
end

struct SortingDirect <: SimulationAlgorithm end

function build_simulator(::SortingDirect, state, model, rates_cache)
  number_jumps = length(model.reactions)
  algorithm = SortingDirectMethod{HasRates}(zeros(number_jumps), 0.0, collect(1:number_jumps), 1)

  return ExactSimulator(algorithm)
end

struct FirstReaction <: SimulationAlgorithm end

function build_simulator(::FirstReaction, state, model, rates_cache)
  number_jumps = length(model.reactions)
  algorithm = FirstReactionMethod(zeros(number_jumps), 0.0)

  return ExactSimulator(algorithm)
end

struct NextReaction <: SimulationAlgorithm end

function build_simulator(::NextReaction, state, model, rates_cache)
  number_jumps = length(model.reactions)
  priority_queue = PQBinaryHeap{Int,Float64,ForwardOrdering}(collect(1:number_jumps), zeros(number_jumps))
  algorithm = NextReactionMethod(zeros(number_jumps), 0.0, priority_queue)

  return ExactSimulator(algorithm)
end
