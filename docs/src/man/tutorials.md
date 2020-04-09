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

```@example neg_autoreg
using DataFrames

DataFrame(ensemble)
```

Each component of the state vector appears as a column, along with trajectory `trial` and timestamp `t` columns.
If the conversion does not produce the expected result, one may be able to force the correct behavior using the `tablefy` function:

```@example neg_autoreg
import BioSimulator: tablefy

DataFrame(tablefy(sample_path)) # DataFrame(sample_path) is currently incorrect
```

## Saving simulation results

Because `SamplePath` and `Ensemble` support iteration, you can save simulation data directly using Julia's I/O interface.
Alternatively, you can use existing packages such as CSV.jl:

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

## Customizing output

It is possible to customize the data stored in a `SamplePath` with the keyword argument `save_function` in [`simulate`](@ref).
This function must accept exactly three arguments:

```
function myfunction(simulator, state, model)
    # your code here
end
```

It has access to data stored in a `simulator` object, the process's `state`, and the underlying `model`.
This customized function should return the data to be stored in a `SamplePath`.
Data can be any object, including custom Julia types.
As an example, the function `save_state` is one of the default options that simply saves a copy of the system state at a given time.
Note the use of annotations to dispatch on the type of `state`

```
# this gets called for well-mixed simulations
save_state(simulator, state::Vector{T}, model) where T <: Int = copy(state)

# this gets called for interacting particle systems
save_state(simulator, state::Lattice, model) = Configuration(state)
```

Another option is `save_rates` which copies the rates for each possible jump

```
save_rates(simulator, state, model) = copy(jump_rates(simulator))
```

BioSimulator provides the following high-level interface for accessing data:

```@docs
cumulative_intensity
jump_rates
next_jump_index
next_jump_time
```

!!! tip

    1. Avoid heavy computations inside a save function. It is usually better to store the data needed and do the work outside the simulation loop. Benchmark different implementations to see what approach works best for your application.

    2. Some data, such as `jump_rates(simulator)` and `state` are stored as numerical arrays. These objects are *mutable* so you should generally return copies rather than the object itself.

    3. The generated `SamplePath` object will support multi-dimensional indexing provided its data is stored as an array. For example, if we set `save_function = save_rates`, then `trajectory[k,j]` will return `rate[j]` at time `trajectory.t[k]`.
