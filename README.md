BioSimulator.jl
===============

BioSimulator is a stochastic simulation package for biochemical reaction networks.

## API Introduction

To illustrate a typical workflow in BioSimulator, we will implement a simple birth-death-immigration process (Kendall's Process).

First, we load BioSimulator:

```jl
using BioSimulator
```

Now we can begin describing our modeling problem. The `Species` type is used to describe individual species or populations that appear in a network. For this example, we will define one `Species` named `Particle` with initial population `5`:

```jl
x1 = Species("Particle", 5, true)
x = [x1]
```

The `Bool` value used in the constructor is used to indicate whether we want to track the population in simulations. The `Species` is then used to create a single entry `Vector`.

Now we define the three reactions that can occur in our network - birth, death, and immigration - using the `Reaction` type:

```jl
r1 = Reaction("Birth", "alpha", [1], [2])
r2 = Reaction("Death", "mu", [1], [0])
r3 = Reaction("Immigration", "nu", [0], [1])
r = [r1, r2, r3]

p = Dict{ASCIIString,Float64}("alpha" => 2.0, "mu" => 1.0, "nu" => 0.5)
```

The first two arguments to the constructor are `String`s identifying the reaction and labeling its reaction rate, respectively. The reaction rate labels must be assigned values in a `Dict` that stores parameters in the network as shown above. The third and fourth arguments are `Vector`s that specify the stoichiometric coefficients of each species as reactants and products, respetively. For example, in the first `Reaction`, one `Particle` occurs as a reactant to produce two `Particle`s. We store the reactions in a `Vector` as before.

Now we can create our `Network`:

```jl
nwk = Network("Kendall's Process", x, r, p)
```

A `Network` can be used to create a `Simulation`:

```jl
sim = Simulation(nwk)
```

A `Simulation` is then passed into a function implementing a stochastic algorithm. BioSimulator currently implements the following algorithms:

* Stochastic Simulation Algorithm (SSA)
* Step Anticipation $\tau$-leaping (SAL)

We can simulate our model using SAL and see the results as follows:

```jl
result1 = sal_explicit(sim, 4.0, tracing=true)
plot_trajectory(result1.traces["Particle"])
```

This will result in an image similar to this:

(TODO)
