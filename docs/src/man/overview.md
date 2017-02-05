# Overview

## Assumptions

Given a system of interacting particles and reactions, we assume:

- the entire particle population is well-mixed (a reaction *will* happen),
- each reaction depends only on the current system state,
- reaction propensities follow the law of mass action.

## Continuous-Time Markov Chains

The underlying mathematical structure in the present modeling approach is a continuous-time Markov chain $\mathbf{X}_{t}$. The random vector $\mathbf{X}_{t}$ has non-negative components $X_{ti}$ that count the number of particles of type $i$ at time $t$. A reaction channel $j$ is characterized by its propensity function $r_{j}(\mathbf{X}_{t})$ that depends on the current vector counts. In a small interval of length $\Delta t$, we expect $r_{j}(\mathbf{x})\Delta t + o(\Delta t)$ reactions of type $j$ to occur. Each reaction $j$ changes the system state by a fixed integer vector $\mathbf{v}^{j}$.

From the Markov Chain perspective, the chain waits an exponential length of time proportional to the cumulative intensity $r_{0}(\mathbf{x}) = \sum_{j=1}^{c} r_{j}(\mathbf{x})$. Once a jump occurs, a reaction channel $j$ fires with probability $r_{j}(\mathbf{x})/r_{0}(\mathbf{x})$.

Let $p_{xy}(t)$ denote the finite-time transition probability of going from state $x$ at time $0$ to state $y$ at time $t$. The probability of transition from $x$ to $y$ during a small interval $\Delta t$ can be approximated to leading order:

\\[ p_{xy}(t + \Delta t) = p_{xy}(t)\left[ 1 - \sum_{j=1}^{c} r_{j}(y) \Delta t \right] + \sum_{j=1}^{c} p_{x,y-v^{j}}(t) r_{j}(y-v^{j}) \Delta t + o(\Delta t). \\]

One can then form a difference quotient and take the limit $\Delta t \to 0$ to arrive at the master equation:

\\[ \frac{\mathrm{d}}{\mathrm{d}t} p_{xy}(t) = \sum_{j=1}^{c} \left[ p_{x,y-v^{j}}(t) r_{j}(y - v^{j}) - p_{xy}(t) r_{j}(y) \right], \\]

which forms a potentially infinite system of coupled ordinary differential equations. The core of `BioSimulator` and stochastic simulation in general is to sample the master equation to arrive at statistically exact results.

## Example Reaction Channels

|      *Name*      |              *Reaction*           | $r(\boldsymbol{x})$ |    $v$    |
| ----------------:| --------------------------------- | ------------------- | --------- |
| Immigration      | $\emptyset \to S_{1}$             | $a_{1}$                  | $[ 1\  0\  0\  0]$ |
| Decay            | $S_{1} \to \emptyset$             | $a_{2}x_{1}$             | $[{-1}\ 0\  0\  0]$ |
| Dimerization     | $S_{1} + S_{1} \to S_{2}$         | $a_{3} \binom{x_{1}}{2}$ | $[-2\  1\  0\  0]$ |
| Isomerization    | $S_{1} \to S_{2}$                 | $a_{4} x_{1}$            | $[{-1}\  1\  0\  0]$ |
| Dissociation     | $S_{2} \to S_{1} + S_{1}$         | $a_{5} x_{2}$            | $[ 2\ {-1}\  0\  0]$ |
| Budding          | $S_{1} \to S_{1} + S_{1}$         | $a_{6} x_{1}$            | $[ 1\  0\  0\  0]$ |
| Replacement      | $S_{1} + S_{2} \to S_{2} + S_{2}$ | $a_{7} x_{1} x_{2}$      | $[{-1}\  1\  0\  0]$ |
| Complex Reaction | $S_{1} + S_{2} \to S_{3} + S_{4}$ | $a_{8} x_{1} x_{2}$      | $[{-1}\  {-1}\  1\  1]$ |
