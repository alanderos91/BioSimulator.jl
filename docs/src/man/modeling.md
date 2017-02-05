# Modeling

```@meta
CurrentModule = BioSimulator
```
## Creating a Model

The `Network` type is the starting point of the modeling process in `BioSimulator`. It represents a collection related `Species` that interact through the rules defined by `Reaction`s.

### Interface

```@docs
  Network
```

```@docs
  Species
```

```@docs
  Reaction
```

## Running Simulations

```@docs
  simulate
```

## Visualizing Results

Simulation output is stored as a 3D array within a `PartialHistory` type where

- the *i* dimension spans the species in the `Network`,
- the *j* dimension spans the epochs sampled, and
- the *k* dimension spans the number of trials.

`BioSimulator` provides quick visualizations of sample and mean trajectories, distributions, and phase plots via the `Plots` package. Please see the [Examples](examples.html) page for specific examples.

## Manipulating and Exporting Data

Manipulating the simulation data directly can be cumbersome. `BioSimulator` uses the `DataFrames` package to extract results as an easy-to-use `DataFrame`. We note that `Plots` supports plotting `DataFrame` objects. One can write simulation results to disk directly, but we recommend using the interface provided by `DataFrames`.

## GUI Generation

`BioSimulator` provides **experimental** automatic graphical user interface generation. One must specify the parameters to expose in the interface and provide default values where necessary.
