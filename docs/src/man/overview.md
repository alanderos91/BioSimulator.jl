# Overview

## Assumptions

Given a system of interacting particles and reactions, we assume:

- the entire particle population is well-mixed (a reaction *will* happen),
- each reaction depends only on the current system state, and
- reaction propensities follow the law of mass action,

## Continuous-Time Markov Chains

The underlying mathematical structure in the present modeling approach is a continuous-time Markov chain $\mathbf{X}_{t}$. The random vector $\mathbf{X}_{t}$ has non-negative components $X_{ti}$ that count the number of particles of type $i$ at time $t$. A reaction channel $j$ is characterized by its propensity function $r_{j}(\mathbf{X}_{t})$ that depends on the current vector counts. Each reaction changes the system state by a fixed integer vector $\mathbf{v}^{j}$.
