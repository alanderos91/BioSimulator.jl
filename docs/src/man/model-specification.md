# Overview

BioSimulator.jl simulates stochastic dynamical systems using Markov chain theory.
Because transition probabilities are usually unavailable and infinitesimal generators are often impractical, we provide a few modeling objects to represent a chain.
Two cases are handled: well-mixed and lattice-based systems.
Here the term "well-mixed" is taken to mean that the underlying space for a chain concerns discrete quantities and/or vectors (non-spatial), whereas lattice-based systems model configurations of interacting particles (spatial).

# Well-Mixed Systems

A `Network` object summarizes interactions between populations.
Each `Species` participates in at least one `Reaction`.

# Lattice/Cellular Automata

Interacting Particle Systems (IPSs) are defined using the `@def_reactions` macro.
The syntax resembles chemical reaction notation but the `+` symbol implies directionality; for example:

```julia
interactions = @def_reactions begin
    X + 0 --> 0 + X, α
    X + Y --> Z + 0, β
    Z --> 0, α
end α β
```

The first interaction reads as "`X` moves to an empty site" (swapping locations).
The second is taken to mean "`X` interacts with an adjacent `Y` to form `Z` at the same location as `X`, with `Y` being consumed in the process".
However, this construction not contain any spatial information.
You can specify a topology using `@enumerate_with_sclass`:

```
model = @enumerate_with_sclass interactions VonNeumann() 2
```

This command "compiles" a 2-D version of the IPS in which interactions are between nearest-neighbors along cardinal directions.
Initial conditions are specified with a `Lattice` object:

```
state = Lattice(coordinates, ...)
```
