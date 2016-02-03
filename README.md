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
m <= Species(:X, 5, istracked=true)
```
The `Network` type overloads the `<=` operator to streamline the model description process. The `Species` type accepts two arguments - an identifier and an initial copy number - as well as an optional argument (`istracked`) indicating whether to track the `Species` through a simulation. This optional argument defaults to `true`.

Now we define the three reactions that can occur in our network - birth, death, and immigration - using the `Reaction` type:

```jl
m <= Reaction(:Birth, :α, r=(:X => 1), p=(:X =>2))
m <= Reaction(:Death, :μ, r=(:X => 1))
m <= Reaction(:Immmigration, :ν, p=(:X => 1))
```

The first argument is an identifier for a reaction, while the second argument assigns an identifier for the rate. The keyword arguments `r` and `p` stand for reactants and products, respectively. A reactant or product list is enclosed in parentheses. These lists consist of (`Symbol`, `Integer`) pairs representing a `Species` and a stoichiometric coefficient.

Next we create the `Parameter`s referenced in our reaction definitions:

```jl
m <= parameter(:α, 2.0, description="birth rate")
m <= parameter(:μ, 1.0, description="death rate")
m <= parameter(:ν, 0.5, description="immigration rate")
```

Each `Parameter` requires an identifier and value. An optional `description` string may be provided to help annotate the model.

Finally, we can simulate the model:

```jl
result = simulate(m, tf=4.0, with=:sal, output=Uniform(), dt=0.1, itr=100_000)
```

This will run 100,000 realizations of our model until simulation time reaches 4 seconds using the Step Anticipation τ-Leaping algorithm (citation). The system's state is recorded every 0.1 seconds. The output returned is a `DataFrame` (from the `DataFrames` package) with columns `Time` and `X`.

We can plot the mean trajectory using the `DataFrames` and `Gadfly` packages:

```jl
using DataFrames
using Gadfly

mean_traj = aggregate(result, :Time, mean)
plot(mean_traj, x=:Time, y=:X_mean, Geom.line)
```

## Supported Algorithms

### Exact Algorithms
* :ssa - Stochastic Simulation Algorithm (Gillespie Algorithm/Direct Method)
* :odm - Optimized Direct Method
* :nrm - Next-Reaction Method
* :frm - First Reaction Method

### Approximate Algorithms
* :sal - Step Anticipation τ-leaping

## Simulation Output

* `Explicit`  - Record network state every time the algorithm steps through, at every iteration.
* `Uniform`   - Record network state at evenly spaced intervals, at every iteration (requires `dt` keyword argument).
* `Histogram` - Record network state at the end of each iteration.
