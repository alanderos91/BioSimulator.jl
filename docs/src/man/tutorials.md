```@setup neg_autoreg
using BioSimulator

pkgdir = BioSimulator |> pathof |> dirname |> dirname
include(joinpath(pkgdir, "test", "test-models", "autoreg.jl"))

network = neg_autoreg();
state, model = BioSimulator.parse_model(network);
```

## Accessing simulation data

A single run generates a `SamplePath` object which is a recursive array provided by [RecursiveArrayTools.jl](https://github.com/JuliaDiffEq/RecursiveArrayTools.jl).

```@example neg_autoreg
sample_path = simulate(state, model, SortingDirect(), tfinal = 100.0, save_points = 0:25:100.0)

sample_path
```

Access the state at the first time point:

```@example neg_autoreg
sample_path[1]
```

Access the third component at the second time point:

```@example neg_autoreg
sample_path[3,2]
```

Access entire history for the third and fourth components:

```@example neg_autoreg
sample_path[[3,4], :]
```

Running multiple simulations generates an `Ensemble`, a collection of `SamplePath` objects; that is, `Vector{SamplePath}`.
These objects also support indexing:

```@example neg_autoreg
# ensemble of 10 trajectories
ensemble = [simulate(state, model, SortingDirect(), tfinal = 100.0, save_points = 0:25:100.0) for _ in 1:10];
```

Retrieve an individual sample path:
```@example neg_autoreg
ensemble[2]
```

Index into a SamplePath:
```@example neg_autoreg
ensemble[2][[3,4], :]
```

## Working with `DataFrames`

Both `SamplePath` and `Ensemble` implement the [IterableTables.jl](https://github.com/queryverse/IterableTables.jl) interface.
This means that they act as *data sources* that can be converted into any supported *data sink*.
For example, you can convert an `Ensemble` into a `DataFrame`:

```
using DataFrames

DataFrame(ensemble)
```

Each component of the state vector appears as a column, along with trajectory `trial` and timestamp `t` columns.
If the conversion does not produce the expected result, one may be able to force the correct behavior using the `tablefy` function:

```
import BioSimulator: tablefy

DataFrame(tablefy(sample_path)) # DataFrame(sample_path) is currently incorrect
```

## Saving simulation results

Because `SamplePath` and `Ensemble` support iteration, you can save simulation data directly using Julia's I/O interface.

> Example here

The easiest approach takes advantage of IterableTables.jl:

```
using CSV

CSV.write(file, result, delim = '\t')
```

## Obtaining summary statistics

Summary statistics for ensembles are supported directly:

```
using Statistics

mean(ensemble)
std(result)
var(result)
```
