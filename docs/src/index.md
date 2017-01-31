# BioSimulator.jl

*A stochastic simulation framework for Julia.*

## Introduction

Many complex systems in biology are analytically intractable, and dynamical predictions based on deterministic models can be grossly misleading. Stochastic simulation algorithms based on continuous-time Markov chains allow researchers to generate accurate time-evolution trajectories, test the sensitivity of models to key parameters, and quantify frequencies of rare events [6, 18, 13].

Situations where stochastic simulation is especially helpful involve:

- rare events such as extinction and mutation,
- key molecules present in small numbers,
- rare reactions with dramatic influence on the dynamics of the system
- population cycles arising from demographic stochasticity.

Examples of such systems include gene expression networks, tumor suppressor pathways, and demographic and ecological systems.

`BioSimulator.jl` aims to provide researchers interested in such phenomena with a fast, reliable, user-friendly, and open-source modeling tool in Julia.

## Installation

`BioSimulator` is not yet registered and must be installed with `Pkg.clone`:

```julia
Pkg.clone("https://alanderos@bitbucket.org/alanderos/biosimulator.jl.git")
```

This package supports Julia `0.5`.

## Table of Contents

##### Overview
##### Modeling
##### Examples
##### Benchmarks
##### Developers
