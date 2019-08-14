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

##### tau-leaping

struct TauLeapingDG2001 <: SimulationAlgorithm end

function build_simulator(::TauLeapingDG2001, state, model, rates_cache)
  number_species = length(state)
  number_jumps = get_number_jumps(model)

  rates = zeros(number_jumps)
  total_rate = zero(eltype(rates))
  leap_formula = DG2001Eq26a(number_species, number_jumps, 0.125)

  # build closure to apply leap updates
  V = extract_net_stoichiometry(model)
  execute_leap! = ApplyLeapUpdate(forward_leap!, V)  # apply a leap update
  reverse_leap! = ApplyLeapUpdate(backward_leap!, V) # reverse a leap update

  # build closure that ensures leaps are valid
  rejection_threshold = 0.75
  proposal = copy(state)
  validate_leap! = RejectionThinning(rejection_threshold, proposal, execute_leap!, reverse_leap!)

  algorithm = PoissonLeapMethod(rates, total_rate, V, leap_formula, validate_leap!)

  return TauLeapSimulator(algorithm, number_jumps, execute_leap!)
end


struct TauLeapingDGLP2003 <: SimulationAlgorithm end

function build_simulator(::TauLeapingDGLP2003, state, model, rates_cache)
  number_species = length(state)
  number_jumps = get_number_jumps(model)

  rates = zeros(number_jumps)
  total_rate = zero(eltype(rates))

  # extract stoichiometry
  V = extract_net_stoichiometry(model)

  # build leap formula
  leap_formula = DGLP2003Eq6(number_species, number_jumps, V, 0.125)

  # build closure to apply leap updates
  execute_leap! = ApplyLeapUpdate(forward_leap!, V)  # apply a leap update
  reverse_leap! = ApplyLeapUpdate(backward_leap!, V) # reverse a leap update

  # build closure that ensures leaps are valid
  rejection_threshold = 0.75
  proposal = copy(state)
  validate_leap! = RejectionThinning(rejection_threshold, proposal, execute_leap!, reverse_leap!)

  algorithm = PoissonLeapMethod(rates, total_rate, V, leap_formula, validate_leap!)

  return TauLeapSimulator(algorithm, number_jumps, execute_leap!)
end

# struct StepAnticipation <: SimulationAlgorithm end

# function build_simulator(::StepAnticipation, state, model, rates_cache)
#   number_species = length(state)
#   number_jumps = get_number_jumps(model)

#   rates = zeros(number_jumps)
#   total_rate = zero(eltype(rates))

#   # extract stoichiometry
#   V = extract_net_stoichiometry(model)

#   # build leap formula

#   # build closure to apply leap updates
#   execute_leap! = ApplyLeapUpdate(forward_leap!, V)  # apply a leap update
#   reverse_leap! = ApplyLeapUpdate(backward_leap!, V) # reverse a leap update

#   # build closure that ensures leaps are valid
#   rejection_threshold = 0.75
#   proposal = copy(state)
#   validate_leap! = RejectionThinning(rejection_threshold, proposal, execute_leap!, reverse_leap!)

#   algorithm = StepAnticipationMethod(rates, total_rate, V, leap_formula, validate_leap!)

#   return TauLeapSimulator(algorithm, number_jumps, execute_leap!)
# end