# Modeling

```@meta
CurrentModule = BioSimulator
```

## Interface

The `Network` type is the starting point of the modeling process in `BioSimulator`. It represents a collection related `Species` that interact through the rules defined by `Reaction`s.

```@docs
  Network
```

A `Species` is simply a named quantity that represents a population.

```@docs
  Species
```

```@docs
  Reaction
```

## Creating a Model

## Running Simulations

## Visualizing Results

Simulation output is stored within a `PartialHistory` as a 3D array where

- the *x* dimension spans the species in the `Network`,
- the *y* dimension spans the epochs sampled, and
- the *z* dimension spans the number of trials.

`BioSimulator` provides convenience methods to quickly visualize sample and mean trajectories, distributions, and phase plots via the `Plots` package.

## Manipulating and Exporting Data

Manipulating the simulation data directly can be cumbersome. `BioSimulator` uses the `DataFrames` package to extract results as an easy-to-use `DataFrame`. We note that `Plots` supports plotting `DataFrame` objects. One can write simulation results to disk directly, but we recommend using the interface provided by `DataFrames`.

## GUI Generation

`BioSimulator` provides **experimental** automatic graphical user interface generation. One simply specifies the fields that should be exposed in the interface while providing defaults for any outstanding fields.
