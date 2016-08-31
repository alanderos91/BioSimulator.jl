"""
```
PartialHistory(t, data, id2ind)
```

Stores the partial output of a simulation job, discretized to specific intervals according to the `sampling_interval` value used when calling `simulate`. Output is stored as a `(t, data)` tuple, where `t` is the time array and data is a 3D array. For example:

- `data[i, j, k]` accesses the value of species at interval `j` along the `k`-th realization.
- `data[:, :, k]` returns the simulation data for the `k`-th realization
- `data[:, j, :]` returns a histogram of every species during the `j`-th interval (i.e. t[j]).
- `data[i, :, :]` returns all the simulation data for the `i`-th species.

For convenience, one may access the simulation data using a special indexing scheme:

```
result = simulate(...) # returns a PartialHistory
result[:X]  # indexing by symbols returns the simulation data for species `X`
result[5.0] # indexing by floats returns simulation data for every species and realization at time `5.0`.
result[100] # indexing by integers returns simulation data for the 100th realization.
```

"""
type PartialHistory
  t      :: LinSpace{Float64}
  data   :: Array{Int, 3}
  id2ind :: Dict{Symbol, Int}
end

function PartialHistory(
    d1      :: Integer, # number of species
    d2      :: Integer, # number of intervals
    d3      :: Integer, # number of realizations
    start   :: Float64, # start time
    stop    :: Float64, # end time
    id2ind  :: Dict{Symbol,Int}
)

  return PartialHistory(linspace(start, stop, d2), zeros(Int, d1, d2, d3), id2ind)
end

function Base.show(io::IO, x::PartialHistory)
  print(io, "[ Partial Simulation History ]\n")
  print(io, "  * species = ", keys(x.id2ind), "\n")
  print(io, "  * no. realizations = ", size(x.data, 3), "\n")
  print(io, "  * start time = ", x.t.start, "\n")
  print(io, "  * end time = ", x.t.stop, " \n")
  print(io, "  * interval stride = ", x.t.stop / x.t.divisor)
end

@inbounds function update!(
    x           :: PartialHistory,
    Xt          :: Vector,
    t           :: AbstractFloat,
    interval    :: Integer,
    realization :: Integer
)

  ttt      = x.t
  data     = x.data
  maxindex = length(ttt)

  while (interval <= maxindex) && (t >= ttt[interval])
    for k in eachindex(Xt)
      data[k, interval, realization] = Xt[k]
    end
    interval = interval + 1
  end

  return interval
end

# indexing with symbols should return a species' history
function Base.getindex(x::PartialHistory, key::Symbol)
  data = x.data
  id2ind = x.id2ind

  i = id2ind[key]

  return reshape(data[i, : , :], size(data, 2), size(data, 3))
end

# indexing with a number should return a particular realization
function Base.getindex(x::PartialHistory, i::Integer)
  data = x.data

  return reshape(data[:, :, i], size(data, 1), size(data, 2))
end

# indexing with a time (float) should return a histogram
function Base.getindex(x::PartialHistory, tn::AbstractFloat)
  data = x.data

  t = x.t
  i = findfirst(t, tn)

  return reshape(data[:, i, :], size(data, 1), size(data, 3))
end
