BioSimulator.jl
===============

BioSimulator is a stochastic simulation package for biochemical reaction networks.

## How To Use

To illustrate a typical workflow in BioSimulator, we will implement a simple birth-death-immigration process (Kendall's Process).

First, we load BioSimulator:

```jl
using BioSimulator
```

and initialize the `Network`:

```jl
m = Network("Kendall's Process")
```
Now we can begin describing our modeling problem. The `Species` type is used to describe distinct populations (e.g. a certain molecule) in the network. For this example, we will define one `Species` named `X` with initial population `5`:

```jl
m <= Species("X", 5)
```
The `Network` type overloads the `<=` operator to streamline the model description process. The `Species` type accepts two arguments - an identifier and an initial copy number.

Next, we define the three reactions that can occur in our network - birth, death, and immigration - using the `Reaction` type:

```jl
m <= Reaction(:Birth, :α, :(X --> X + X))
m <= Reaction(:Death, :μ, :(X --> 0)
m <= Reaction(:Immmigration, :ν, :(0 --> X))
```

The first argument is an identifier for a reaction, while the second argument assigns an identifier for the rate. The last argument is a `Expr` object that mimics a chemical formula. Reactants appear on the left, while products appear on the right. A `0` is used to denote the absence of a reactant or product. Additionally, one may specify stoichiometric coefficients e.g. `:(X --> 2*X)` or `:(X --> 2X)`.

Lastly, we assign values to the `Parameter`s referenced in our reaction definitions:

```jl
m <= parameter(:α, 2.0, description="birth rate")
m <= parameter(:μ, 1.0, description="death rate")
m <= parameter(:ν, 0.5, description="immigration rate")
```

Each `Parameter` requires an identifier and value. An optional `label` may be provided to help annotate the model.

Finally, we can simulate the model:

```jl
result = simulate(m, time=4.0, method=:SAL, output=:explicit, sampling_interval, realizations=100_000)
```

This will run 100,000 realizations of our model until simulation time reaches 4 units of time using the Step Anticipation τ-Leaping algorithm (Sehl et al). The system's state is recorded every 0.1 units of time. The output returned is a `SimulationOutput` object. It stores the time series data as a `DataFrame`. Additionally, a `Dict` stores algorithm statistics and other simulation metadata of interest.

## Supported Algorithms

### Exact Algorithms
* :SSA - Stochastic Simulation Algorithm (Gillespie Algorithm/Direct Method)
* :ODM - Optimized Direct Method
* :NRM - Next-Reaction Method
* :FRM - First Reaction Method

### Approximate Algorithms
* :SAL - Step Anticipation τ-leaping

## Simulation Output

* `:explicit`  - Record network state every time the algorithm steps through, at every iteration.
* `:fixed`   - Record network state at evenly spaced intervals, at every iteration (requires `sampling_interval` keyword argument).
