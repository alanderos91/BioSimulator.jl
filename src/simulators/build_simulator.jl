# dirty hacks
get_number_jumps(model::ReactionSystem) = length(model.reactions)
get_number_jumps(model::InteractingParticleSystem) = length(model.reactions)
get_number_jumps(model::Vector{T}) where T = length(model)
get_number_species(state::Lattice) = length(state.types)
get_number_species(state::Vector{T}) where T = length(state)

abstract type SimulationAlgorithm end

"""
```
Direct()
```

Gillespie's Direct Method.

#### References

- Gillespie, D.T. (1976) A general method for numerically simulating the stochastic time evolution of coupled chemical reactions. *Journal of Computational Physics*. [https://doi.org/10.1016/0021-9991(76)90041-3](https://doi.org/10.1016/0021-9991(76)90041-3)

- Gillespie, D.T. (1977) Exact stochastic simulation of coupled chemical reactions. *The Journal of Physical Chemistry*. [https://doi.org/10.1021/j100540a008](https://doi.org/10.1021/j100540a008)
"""
struct Direct <: SimulationAlgorithm end

function build_simulator(::Direct, state, model, rates_cache)
  number_jumps = get_number_jumps(model)
  algorithm = DirectMethod{rates_cache}(zeros(number_jumps), 0.0)

  return ExactSimulator(algorithm)
end

"""
```
EnhancedDirect()
```

A Direct method that uses a dependency graph.

#### References

- Marchetti L., Priami C., Thanh V.H. (2017) Simulation Algorithms for Computational Systems Biology. *Texts in Theoretical Computer Science*.
"""
struct EnhancedDirect <: SimulationAlgorithm end

function build_simulator(::EnhancedDirect, state, model, rates_cache)
  number_jumps = get_number_jumps(model)
  algorithm = EnhancedDirectMethod{HasRates}(zeros(number_jumps), 0.0)

  return ExactSimulator(algorithm)
end

"""
```
SortingDirect()
```

The sorting direct method reorders reactions to identify frequently occurring events.
This effectively reduces the search depth in selecting the next reaction.

#### References

- McCollum, J.M., Peterson, G.D., Cox, C.D., Simpson, M.L., Samatova, N.F. (2006). The sorting direct method for stochastic simulation of biochemical systems with varying reaction execution behavior. *Computational Biology and Chemistry*. [https://doi.org/10.1016/j.compbiolchem.2005.10.007](https://doi.org/10.1016/j.compbiolchem.2005.10.007).
"""
struct SortingDirect <: SimulationAlgorithm end

function build_simulator(::SortingDirect, state, model, rates_cache)
  number_jumps = get_number_jumps(model)
  algorithm = SortingDirectMethod{HasRates}(zeros(number_jumps), 0.0, collect(1:number_jumps), 1)

  return ExactSimulator(algorithm)
end

"""
```
FirstReaction()
```

A first reaction method that samples reaction events and times jointly.

#### References

- Gillespie, D.T. (1976) A general method for numerically simulating the stochastic time evolution of coupled chemical reactions. *Journal of Computational Physics*. [https://doi.org/10.1016/0021-9991(76)90041-3](https://doi.org/10.1016/0021-9991(76)90041-3)
"""
struct FirstReaction <: SimulationAlgorithm end

function build_simulator(::FirstReaction, state, model, rates_cache)
  number_jumps = get_number_jumps(model)
  algorithm = FirstReactionMethod(zeros(number_jumps), 0.0)

  return ExactSimulator(algorithm)
end

"""
```
NextReaction()
```

Gibson and Bruck's next reaction method.
Uses a priority queue to select the next reaction event and time.

#### References

- Gibson, M.A., Bruck, J. (1999) Efficient exact stochastic simulation of chemical systems with many species and many channels. *Journal of Physical Chemistry*. [https://doi.org/10.1021/jp993732q](https://doi.org/10.1021/jp993732q)
"""
struct NextReaction <: SimulationAlgorithm end

function build_simulator(::NextReaction, state, model, rates_cache)
  number_jumps = get_number_jumps(model)
  priority_queue = PQBinaryHeap{Int,Float64,ForwardOrdering}(collect(1:number_jumps), zeros(number_jumps))
  algorithm = NextReactionMethod(zeros(number_jumps), 0.0, priority_queue)

  return ExactSimulator(algorithm)
end

"""
```
RejectionSSA()
```

A rejection-based algorithm that operates on "virtual states" to accelerate simulation.

WARNING: This *implementation* is a work in progress.

#### References

- Thahn, V.H., Priami, C., Zunino, R. (2014) Efficient rejection-based simulation of biochemical reactions with stochastic noise and delays. *Journal of Chemical Physics*. [https://doi.org/10.1063/1.4896985](https://doi.org/10.1063/1.4896985)
"""
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

"""
```
TauLeapingDG2001()
```

A pure tau-leaping method using the Equation (26a) to select tau.

#### Optional Arguments

Work in progress.

#### References

- Gillespie, D.T. (2001) Approximate accelerated stochastic simulation of chemically reacting systems. *Journal of Chemical Physics*. [https://doi.org/10.1063/1.1378322](https://doi.org/10.1063/1.1378322)
"""
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

"""
```
TauLeapingDGLP2003()
```

A pure tau-leaping method using the Equation (6) to select tau.

#### Optional Arguments

Work in progress.

#### References

- Gillespie, D.T., Petzold, L.R. (2003) Improved leap-size selection for accelerated stochastic simulation. *Journal of Chemical Physics*. [https://doi.org/10.1063/1.1613254](https://doi.org/10.1063/1.1613254)
"""
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


"""
```
StepAnticipation()
```

A pure tau-leaping method using the Equation (15) to select tau.

#### Optional Arguments

Work in progress.

#### References

- Sehl, M.E., Alekseyenko, A.L., Lange, K.L. (2009) Accurate stochastic simulation via the step anticipation tau-leaping (SAL) algorithm. *Journal of Computational Biology*. [https://dx.doi.org/10.1089/cmb.2008.0249](https://dx.doi.org/10.1089/cmb.2008.0249)
"""
struct StepAnticipation <: SimulationAlgorithm end

function build_simulator(::StepAnticipation, state, model, rates_cache)
  number_species = length(state)
  number_jumps = get_number_jumps(model)

  rates = zeros(number_jumps)
  total_rate = zero(eltype(rates))

  # extract stoichiometry
  U = extract_coefficients(model)
  V = extract_net_stoichiometry(model)

  # allocate memory for derivatives
  dxdt = zeros(number_species)
  drdt = zeros(number_jumps)

  # build leap formula
  k = model.rxn_rates # note: this needs to change in the future!
  leap_formula = GenericLeapFormula(k, U, V, dxdt, drdt, 0.125)

  # build closure to apply leap updates
  execute_leap! = ApplyLeapUpdate(forward_leap!, V)  # apply a leap update
  reverse_leap! = ApplyLeapUpdate(backward_leap!, V) # reverse a leap update

  # build closure that ensures leaps are valid
  rejection_threshold = 0.75
  proposal = copy(state)
  validate_leap! = RejectionThinning(rejection_threshold, proposal, execute_leap!, reverse_leap!)

  algorithm = StepAnticipationMethod(rates, total_rate, leap_formula, validate_leap!, U, V, dxdt, drdt)

  return TauLeapSimulator(algorithm, number_jumps, execute_leap!)
end

"""
```
HybridSAL()
```

Same as `StepAnticipation()`, but defaults to `Direct()` depending on the cumulative intensity.

#### Optional Arguments

Work in progress.

#### References

- Sehl, M.E., Alekseyenko, A.L., Lange, K.L. (2009) Accurate stochastic simulation via the step anticipation tau-leaping (SAL) algorithm. *Journal of Computational Biology*. [https://dx.doi.org/10.1089/cmb.2008.0249](https://dx.doi.org/10.1089/cmb.2008.0249)
"""
struct HybridSAL <: SimulationAlgorithm end

function build_simulator(::HybridSAL, state, model, rates_cache)
  # shared information
  number_species = length(state)
  number_jumps = get_number_jumps(model)

  rates = zeros(number_jumps)
  total_rate = zero(eltype(rates))

  # exact simulator...
  algorithm = DirectMethod{rates_cache}(rates, 0.0)
  exact = ExactSimulator(algorithm)

  # tau-leaping simulator...

  # extract stoichiometry
  U = extract_coefficients(model)
  V = extract_net_stoichiometry(model)

  # allocate memory for derivatives
  dxdt = zeros(number_species)
  drdt = zeros(number_jumps)

  # build leap formula
  k = model.rxn_rates # note: this needs to change in the future!
  leap_formula = GenericLeapFormula(k, U, V, dxdt, drdt, 0.125)

  # build closure to apply leap updates
  execute_leap! = ApplyLeapUpdate(forward_leap!, V)  # apply a leap update
  reverse_leap! = ApplyLeapUpdate(backward_leap!, V) # reverse a leap update

  # build closure that ensures leaps are valid
  rejection_threshold = 0.75
  proposal = copy(state)
  validate_leap! = RejectionThinning(rejection_threshold, proposal, execute_leap!, reverse_leap!)

  algorithm = StepAnticipationMethod(rates, total_rate, leap_formula, validate_leap!, U, V, dxdt, drdt)

  tauleap = TauLeapSimulator(algorithm, number_jumps, execute_leap!)

  return HybridTauLeapSimulator(exact, tauleap)
end
