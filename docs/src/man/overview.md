# Overview

```@meta
CurrentModule = BioSimulator
```

## Creating a Model

The `Network` type is the starting point of the modeling process in BioSimulator.jl.
A `Network` is a system of coupled `Species` and `Reaction` channels.

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

### Parallelization

One can start Julia with the `--procs=N` option to enable parallelization. Here `N` is the number of threads available.
In this case, the `simulate` routine will automatically delegate individual simulations to each available thread until the total number of `trials` is satisfied.

## Simulation results

The result of `simulate` is a `SimulationSummary` that stores information about a model, the algorithm used in a simulation, and simulation results.
We summarize the fields of a `result` of type `SimulationSummary` below:

* `result.model`: The `Network` that was simulated.
* `result.algorithm_name`: The name of the algorithm used in the simulation.
* `result.algorithm_params`: The settings for the algorithm, including the epochs, time span, and trials used.
* `result.algorithm_stats`: Any statistics for the algorithm used. Direct methods simply record the number of steps taken over all trials. Tau-leaping methods also record the number of negative excursions and the number of times the recovery mechanism is used.
* `result.simulation_data`: The results of the simulation. It is a vector whose elements are `SamplePaths` (using the `Val{:full}` option) or `RegularPaths` (using the `Val{:fixed}` option).
* `result.id2index`: A dictionary that maps the name of a species to an index into the state vector $\mathbf{X}_{t}$.

### SamplePath output

Using the `Val{:full}` option for `simulate` records the state vector of the process after every event.
In this case, the `simulation_data` field is a vector of `SamplePath`s.
Each entry of the vector is an individual realization of the process.
A `SamplePath` is a variable-length collection of records at various time points.
Given a `SamplePath`, say `xw`, one can access the time points by `xw.tdata` and the states by `xw.xdata`.
Specifically, the `i`-th record corresponds to a state $\mathbf{X}_{t_{i}}$ at time $t_{i}$.

The advantage of `SamplePath` output and the `Val{:full}` option is that one can recover the distribution of each species at different points in time.

### RegularPath output

Using the `Val{:fixed}` option for `simulate` records the state vector at pre-defined intervals determined by the `epochs` option.
Given that `epochs = n`, the time span of a simulation is broken up into `n` intervals of equal length.
After each step, the simulation checks to see if it has stepped over into the next epoch.
If it has, then it records the last state as the value of the stochastic process at the *previous* epoch.
This option produces a vector of `RegularPath`s, with each entry corresponding to an independent realization of the process.
Each `RegularPath` has the fields `tdata` and `xdata` representing time points and states, respectively.
Importantly, `xdata` is a pre-allocated matrix with column `i` storing the value $\mathbf{X}_{t_{i}}$.

One advantage of `RegularPath` output and the `Val{:fixed}` option is a slight performance boost attributed to minimal dynamic memory allocation.
This advantage may disappear if one uses far more `epochs` than the expected number of events in a given time span.
However, care must be taken in choosing the number of epochs.
Using too few epochs may grossly simplify the observed qualitative behavior in visualizations.
Moreover, only the last save point may be used to reliably estimate a probability distribution.
Another advantage of this option is that `RegularPath`s facilitate computing summary statistics, including mean trajectories.

## Plotting

BioSimulator.jl provides plot recipes that operates on a `SimulationSummary`.
The basic interface is as follows:

```julia
plot(result::SimulationSummary{T}; plot_type = :trajectory,
    trial = nothing, species = nothing, epochs = nothing) where T
```

There are three different options for `plot_type`: `:trajectory`, `:meantrajectory`, and `:histogram`.
Each of these options may be modulated by the `trial`, `species`, and `epochs` options.

### Trajectory

The `plot_type = :trajectory` option produces the time evolution of a single realization of the underlying stochastic process.
The `trial` option can be used to select different realizations.
The `species` option can be used to specify a vector of species to plot in a figure.

### Mean Trajectory

The `plot_type = :meantrajectory` option produces the *averaged* time evolution of the underlying stochastic process.
The expected value of a particular species is plotted as a point at evenly spaced intervals.
Error bars denote one standard deviation away from the mean.
Note that this option is fully compatible with the `Val{:full}` option by discretizing the time span.
This is achieved by passing a value for the `epochs` option.
The `species` option can be used to specify a vector of species to plot in a figure.

### Histogram

The `plot_type = :histogram` option produces an estimate of the distribution of the process at the final time point in the simulation.
The `species` option can be used to specify a vector of species to plot in a figure.
