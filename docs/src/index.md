# BioSimulator.jl

*A stochastic simulation framework for Julia.*

## Introduction

## Installation

`BioSimulator` is not yet registered and must be installed with `Pkg.clone`:

```julia
Pkg.clone("https://alanderos@bitbucket.org/alanderos/biosimulator.jl.git")
```

This package supports Julia `0.5`.

## Table of Contents

##### Stochastic Simulation Overview
```@contents
Pages = [
  "man/modeling_assumptions.md",
  "man/algorithms.md"
]
Depth = 1
```
##### Modeling
```@contents
Pages = [
  "man/interface.md",
  "man/creating_a_model.md",
  "man/running_simulations.md",
  "man/visualizing_results.md",
  "man/gui_generations.md"
]
Depth = 1
```
##### Examples
```@contents
Pages = [
  "man/birth-death-immigration_process.md",
  "man/michaelis-menten_enzyme_kinetics.md",
  "man/prokaryotic_auto-regulatory_gene_network.md"
]
Depth = 1
```
##### Benchmarks
##### Developers
