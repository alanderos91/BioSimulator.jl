## Basics

```@setup neg_autoreg
ENV["GKSwstype"] = "nul"

using BioSimulator, Plots

gr(fmt = :png, dpi = 200)

pkgdir = BioSimulator |> pathof |> dirname |> dirname
include(joinpath(pkgdir, "test", "test-models", "autoreg.jl"))

network = neg_autoreg();
state, model = BioSimulator.parse_model(network);
result = [simulate(state, model, SortingDirect(), tfinal = 500.0, save_points = 0:25:500.0) for _ in 1:100];
```

The plotting defaults provided by BioSimulator.jl require the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) package. You can install it with

```
Pkg.add("Plots")
```

BioSimulator.jl does not load the Plots.jl package by default.
Any time you need plotting functionality, simply load the package:

```
# if BioSimulator is already loaded
using Plots

# if you're just starting
using BioSimulator, Plots
```

Note that Plots.jl is independent of BioSimulator.jl.
Please consult the [Plots.jl documentation](http://docs.juliaplots.org/latest/) for additional details.

---

BioSimulator.jl provides a few plotting recipes to quickly summarize simulation results.
The type of plot is controlled by a `summary` keyword.
The available options include:

- `summary = :trajectory`: An individual `SamplePath`.
- `summary = :mean`: Averaged trajectories based on an `Ensemble`.
- `summary = :histogram`: A histogram based on an `Ensemble`.
- `summary = :phase`: A phase portrait of the stochastic system. The default puts the first component on $x$-axis and the second component on the $y-$axis.

See the following examples on how to use the default recipes.

## Sample trajectories

```@example neg_autoreg
# `result` is an `Ensemble` of 100 trajectories
plot(result[1], summary = :trajectory,
    xlabel = "time", ylabel = "count",
    label = ["gene" "P2_gene" "mRNA" "P" "P2"])
```

## Mean trajectories

```@example neg_autoreg
# select data for protein `P` and dimer `P2` with their indices
plot(result, summary = :mean,
    xlabel = "time", ylabel = "average",
    vars = [4,5], label = ["P" "P2"])
```

## Distributions

```@example neg_autoreg
# 3 figures on top row, 2 on the bottom
mylayout = @layout [a{0.33w} a{0.33w} a{0.33w}
                        a{0.5w} a{0.5w}]

plot(result, summary = :histogram,
    ylabel = "probability",
    timepoint = 500.0, normalize = :probability, bins = :sturges,
    layout = mylayout,
    label = ["gene" "P2_gene" "mRNA" "P" "P2"])
```

## Phase plots

For a `SamplePath`:

```@example neg_autoreg
plot(result[1], summary = :phase,
    vars = (4,5), xlabel = "P", ylabel = "P2", colorbar = true)
```

For an `Ensemble`:

```@example neg_autoreg
plot(result, summary = :phase,
    vars = (4,5), xlabel = "P", ylabel = "P2", colorbar = true)
```

## Configurations (Lattice-based models)

You can visualize individual configurations in 2D or 3D:

```
plot(result[10])
```

A `SamplePath` is visualized by plotting the initial and final configuration:

```
plot(result)
```

## General tips

- The components of a `SamplePath` follow the same order in the `Network` definition. For example, if we define `X`, `Y`, and `Z` in our model, then their indices are `1`, `2`, and `3`, respectively.

- By default, *every* variable will appear in a plot (except for phase portraits). You can select a few variables using the `var` keyword argument.
