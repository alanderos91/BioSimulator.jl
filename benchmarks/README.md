# Benchmarks

This folder contains the benchmark code used in **BioSimulator.jl: Stochastic simulation in Julia**.
It compares BioSimulator.jl, [StochPy](http://stochpy.sourceforge.net/), [StochKit2](https://sourceforge.net/projects/stochkit/), and [Gillespie.jl](https://github.com/sdwfrost/Gillespie.jl) across a few examples. These benchmarks are not exhaustive.

# What are the benchmarks?

These scripts measure the performance of Gillespie's Direct method (SSA) across each software tool.
In addition, the BioSimulator.jl benchmarks include additional results for the algorithms currently implemented (SSA, FRM, NRM, OTL, and SAL) using full simulation output, fixed-interval output, and multi-threading.
There are five basic models included in this suite (see [model definitions](#Model-definitons)).

The *task* that is benchmarked is

```
run-ssa model t_final n_saves n_trial
```

where `run-ssa` is the command that runs SSA for a given software package.
Here `t_final` specifies a time span for an individual realization of the stochastic process, `n_saves` is the number of time points saved for a single realization, and `n_trial` is the number of realizations to simulate.
The time to complete this task is sampled several times to obtain reasonable estimates of summary statistics.
See the section on [benchmark parameters](#Benchmark-parameters) for additional details.

**The benchmark task we have defined reflects how quickly a given tool generates multiple time trajectories of a given stochastic process using SSA.**
"Good performance" depends on how well the following interact:

- internal model representation/construction
- algorithm implementation
- simulation engine(s)
- individual output generation/storage
- ensemble output generation/storage

With the exception of `StochPy`, performance is measured using [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl).

# Minimal requirements

| Software        | Version | Dependencies               |
|-----------------|---------|----------------------------|
| BioSimulator.jl | 0.4     | Julia v0.6.*               |
| StochPy         | 2.3     | Python 2.6+ or 3.4+, NumPy |
| StochKit2       | 2.0.13  |                            |
| Gillespie.jl    | 0.0.2   |                            |

**Note**: The StochPy scripts use `python3` by default.

# Running the benchmarks

## Setup

The main script is contained in `examples.sh`.
It requires `julia`, `python3`, and `ssa` (for StochKit2) to be visible in your `PATH`.

## Individual scripts

BioSimulator.jl, StochPy, and Gillespie.jl each have an associated script for running a benchmark from the command line.
Running BioSimulator.jl with the multithreading option requires setting the `JULIA_NUM_THREADS` environment variable (see `examples.sh` for an example).
It uses the following scripts to run the suite:

#### BioSimulator.jl

Usage: `julia benchmarks_biosimulator.jl [model1] [model2] ...`

#### StochPy

Usage: `python3 stochpy_bench.py [model] [t_final] [n_saves] [seed] [n_trial]`

#### StochKit

Usage: `julia benchmarks_stochkit.jl [model1] [model2] ...`

#### Gilespie.jl

Usage: `julia benchmarks_gillespie.jl [model1] [model2] ...`

## Model definitions

Directories:
- `autoreg`
- `dimer-decay`
- `kendall`
- `mmek`
- `yeast`

Files (within each model directory):

- `biosimulator.jl`
- `gillespie.jl`
- `stochkit2.xml`
- `stochpy.psc`

## Benchmark parameters

With the exception of the StochPy script, all the benchmark scripts source benchmark parameters from `benchmark_parameters.jl`.
This file simply stores the arguments passed to each tool's simulation command.
Parameters are stored in a `ParamSet` which is a 5-tuple with the following structure:

- `t_final`: Each simulation command realizes the simulated stochastic process up to this value.
- `n_saves`: The number of time points to save for each realization. Used only for fixed-interval output.
- `n_trial`: The number of realizations to simulate.
- `n_sample`: The number of times to sample a simulation *task*, as defined above ([What are the benchmarks?](What-are-the-benchmarks?))
- `t_limit`: A time limit for individual benchmarks. A benchmark task will terminate if this limit is exceeded, even if the number of samples is less than `n_sample`.

In addition, a `SEED` value is used to seed each benchmark task.
The `MODELS` variable stores the name of each model included in this benchmark suite.
The `parameters` dictionary is used to associate `MODELS[i]` to a `ParamSet`.

## Results

Directories: `biosimulator-serial`, `biosimulator-parallel`, `gillespie`, `stochkit`, `stochpy`

The output files are structured as follows:

```
----- <model-name-here> -----
trials: <value>       # no. times the task was measured
  mean: <estimate> ms # mean based on <value>
   std: <estimate> ms # standard deviation based on <value>
median: <estimate> ms # median based on <value>
    Q1: <estimate> ms # first quartile based on <value>
    Q3: <estimate> ms # third quartile based on <value>
```

# Important caveats for consideration

The tools selected here vary in output saving strategies, algorithm selection, and model representation.
These notes should inform the interpretation of the benchmark results.
Bold type is used to emphasize a point that is likely to have a large impact on performance.
Please feel free to file an issue if an important point is missing. 

### BioSimulator.jl

- Uses fixed-interval output.
- Generates a dependency graph even if it is not used.

### StochPy

- Saves state after each simulation step. Fixed-interval output is available but it uses StochKit2 internally.
- Saves additional information about reaction propensities.

### StochKit2

- Uses fixed-interval output.
- Uses parallelism by default.
- Selects a variant of SSA based on model information. In our benchmarks, `odm_ssa_small` is used.

### Gillespie.jl

- Saves state after each simulation step.
- Does not generate trajectories in an ensemble like the other tools. **We implement a simple wrapper that uses a `for` loop to mimic this feature.** The output of each run is effectively ignored.