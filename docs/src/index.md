# BioSimulator.jl

*A stochastic simulation framework for Julia.*

## Introduction

Many complex systems in biology are analytically intractable, and dynamical predictions based on deterministic models can be grossly misleading.
Stochastic simulation algorithms based on continuous-time Markov chains allow researchers to generate accurate time-evolution trajectories, test the sensitivity of models to key parameters, and quantify frequencies of rare events.

Situations where stochastic simulation is especially helpful involve:

- rare events such as extinction and mutation,
- key molecules present in small numbers,
- rare reactions with dramatic influence on the dynamics of the system
- population cycles arising from demographic stochasticity.

Examples of such systems include gene expression networks, tumor suppressor pathways, and demographic and ecological systems.

BioSimulator.jl aims to provide researchers interested in such phenomena with a fast, reliable, user-friendly, and open-source modeling tool in Julia.

## Installation

BioSimulator.jl must be installed with `Pkg.clone` in the Julia REPL:

```julia
Pkg.clone("https://github.com/alanderos91/biosimulator.jl.git", "BioSimulator")
```

You can start using BioSimulator.jl in scripts or the REPL with the command:

```julia
using BioSimulator
```

### Additional tools

#### DataFrames

The `SimulationSummary` returned by `simulate` can be converted into a `DataFrame` using the [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl) package.
The conversion is straightforward after installing the package with the command `Pkg.add("DataFrames")`:

```julia
# make sure the package is loaded
using DataFrames

# simulate a model and save the results
result = simulate(model...)

# returns a DataFrame
DataFrame(result)
```

The first two columns indicate the time point and trial of a record (row).
The remaining columns represent the species counts and are labeled using the given names.

#### Plotting

The plotting defaults provided by BioSimulator.jl require the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) package. You can install it with

```julia
Pkg.add("Plots")
```

BioSimulator.jl does not load the Plots.jl package by default.
Any time you need plotting functionality, simply load the package:

```julia
# if BioSimulator is already loaded
using Plots

# if you're just starting
using BioSimulator, Plots
```

Note that Plots.jl is independent of BioSimulator.jl and can be used without BioSimulator.jl.
Please consult the Plots.jl documentation for additional details.

#### Petri nets

BioSimulator.jl should install the [TikzGraphs.jl](https://github.com/sisl/TikzGraphs.jl) package by default.
You can try generating a Petri net in a Jupyter notebook for a `model` using `visualize(model)`.
If you want to generate the figure in a script and save it to a file, install the [TikzPictures.jl](https://github.com/sisl/TikzPictures.jl) package.
You can then save a figure using:

```
import TikzPictures: save

figure = visualize(model)
save(PDF(filename), figure)
```

## Table of Contents

```@contents
  pages = [
    "Home"       => "index.md",
    "Algorithms" => "man/algorithms.md"
  ]
```
