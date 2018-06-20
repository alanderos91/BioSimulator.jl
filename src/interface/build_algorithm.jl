"""
```
Direct()
```

Gillespie's Direct Method. Simulates a system of coupled reactions and species by computing the time to the next reaction and searching on the CMF.

#### References

Daniel T Gillespie. A general method for numerically simulating the stochastic time evolution of coupled chemical reactions. *Journal of Computational Physics*, 1976. [https://doi.org/10.1016/0021-9991(76)90041-3](https://doi.org/10.1016/0021-9991(76)90041-3)

Daniel T. Gillespie. Exact stochastic simulation of coupled chemical reactions. *The Journal of Physical Chemistry*, 1977. [https://doi.org/10.1021/j100540a008](https://doi.org/10.1021/j100540a008)
"""
struct Direct end

function build_algorithm(::Direct, end_time, track_stats; kwargs...)
  SSA(end_time, track_stats)
end

"""
```
FirstReaction()
```

Gillespie's First Reaction Method.

Statistically equivalent to `SSA`, but more computationally expensive: it computes the time to the next reaction as the minimum waiting time relative to the next firing times of each reaction.

#### References

Daniel T Gillespie. A general method for numerically simulating the stochastic time evolution of coupled chemical reactions. *Journal of Computational Physics*, 1976. [https://doi.org/10.1016/0021-9991(76)90041-3](https://doi.org/10.1016/0021-9991(76)90041-3)
"""
struct FirstReaction end

function build_algorithm(::FirstReaction, end_time, track_stats; kwargs...)
  FRM(end_time, track_stats)
end

"""
```
NextReaction()
```

Gibson and Bruck's Next Reaction Method. It provides better computational efficiency on networks with loosely connected reactions.

#### References

Michael A. Gibson and Jehoshua Bruck. Efficient exact stochastic simulation of chemical systems with many species and many channels. *Journal of Physical Chemistry*, 1999. [https://doi.org/10.1021/jp993732q](https://doi.org/10.1021/jp993732q)
"""
struct NextReaction end

function build_algorithm(::NextReaction, end_time, track_stats; kwargs...)
  NRM(end_time, track_stats)
end

"""
```
OptimizedDirect()
```

Optimized Direct Method. Similar to `SSA`, with the added benefit of sorting reactions according to their propensities over time. This improves the search on the CMF. Sorting is done with pre-simulation.

#### References

Yang Cao, Hong Li, and Linda Petzold. Efficient formulation of the stochastic simulation algorithm for chemically reacting systems. *Journal of Chemical Physics*, 2004. [https://doi.org/10.1063/1.1778376](https://doi.org/10.1063/1.1778376)
"""
struct OptimizedDirect end

function build_algorithm(::OptimizedDirect, end_time, track_stats; kwargs...)
  ODM(end_time, track_stats)
end

"""
```
TauLeaping()
```

Poisson tau-leaping. This version implements the description by Gillespie and Petzold (2003) and provides safeguards against negative populations.

#### Optional Arguments

- `epsilon = 0.03`: A parameter controlling the size of leaps. Higher values allow for larger leaps, but may compromise the accuracy of results.
- `delta = 2.0`: A parameter used to switch between a Gillespie update and tau-leaping. If the next tau leap is less than some multiple `delta` of the expected time for one event to occur, then `StepAnticipation` carries out a Gillespie step to avoid negative species populations and save on efficiency.
- `beta = 0.75`: A parameter between 0 and 1 used to contract a tau leap in the event of negative populations. Aggresive contraction (lower values) will bias sample paths.

#### References

Daniel T. Gillespie and Linda R. Petzold. Improved leap-size selection for accelerated stochastic simulation. *Journal of Chemical Physics*, 2003. [https://doi.org/10.1063/1.1613254](https://doi.org/10.1063/1.1613254)
"""
struct TauLeaping end

function build_algorithm(::TauLeaping, end_time, track_stats;
    epsilon::Float64 = 0.03, delta::Float64 = 2.0, beta::Float64 = 0.75, kwargs...)
  OTL(end_time, ϵ, δ, β, track_stats)
end

"""
```
StepAnticipation()
```

Step Anticipation tau-Leaping. An algorithm that improves upon the accuracy of Poisson tau-leaping techniques by incorporating a first-order Taylor expansion. It provides safeguards against negative populations.

#### Optional Arguments

- `epsilon = 0.03`: A parameter controlling the size of leaps. Higher values allow for larger leaps, but may compromise the accuracy of results.
- `delta = 2.0`: A parameter used to switch between a Gillespie update and tau-leaping. If the next tau leap is less than some multiple `delta` of the expected time for one event to occur, then `StepAnticipation` carries out a Gillespie step to avoid negative species populations and save on efficiency.
- `beta = 0.75`: A parameter between 0 and 1 used to contract a tau leap in the event of negative populations. Aggresive contraction (lower values) will bias sample paths.


#### References

Mary E. Sehl, Alexander L. Alekseyenko, and Kenneth L. Lange. Accurate stochastic simulation via the step anticipation tau-leaping (SAL) algorithm. *Journal of Computational Biology*, 2009. [https://dx.doi.org/10.1089/cmb.2008.0249](https://dx.doi.org/10.1089/cmb.2008.0249)

"""
struct StepAnticipation end

function build_algorithm(::StepAnticipation, end_time, track_stats;
  epsilon::Float64 = 0.03, delta::Float64 = 2.0, beta::Float64 = 0.75, kwargs...)
  SAL(end_time, ϵ, δ, β, track_stats)
end
