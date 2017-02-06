# BioSimulator.jl

*A stochastic simulation framework for Julia.*

## Introduction

Many complex systems in biology are analytically intractable, and dynamical predictions based on deterministic models can be grossly misleading. Stochastic simulation algorithms based on continuous-time Markov chains allow researchers to generate accurate time-evolution trajectories, test the sensitivity of models to key parameters, and quantify frequencies of rare events [[6](man/references.html#6), [11](man/references.html#11), [15](man/references.html#15)].

Situations where stochastic simulation is especially helpful involve:

- rare events such as extinction and mutation,
- key molecules present in small numbers,
- rare reactions with dramatic influence on the dynamics of the system
- population cycles arising from demographic stochasticity.

Examples of such systems include gene expression networks, tumor suppressor pathways, and demographic and ecological systems.

`BioSimulator.jl` aims to provide researchers interested in such phenomena with a fast, reliable, user-friendly, and open-source modeling tool in Julia.

This package supports Julia `0.5`.

## Table of Contents

```@contents
  pages = [
    "Home"       => "index.md",
    "Overview"   => "man/overview.md",
    "Modeling"   => "man/modeling.md",
    "Algorithms" => "man/algorithms.md",
    "Examples"   => "man/examples.md",
    "Benchmarks" => "man/benchmarks.md",
    "Developers" => "man/developers.md",
    "References" => "man/references.md"
  ]
```
