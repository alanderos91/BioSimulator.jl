struct SimData{T}
  id2ind  :: Dict{Symbol,Int}
  t_index :: StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
  data    :: T
  stats   :: Dict{Symbol,Int}
end

get_data(output :: SimData) = output.t_index, output.data

@inline function Base.getindex(output :: SimData, key :: Symbol)
  index = output.id2ind[key]
  return output.data[index, :, :]
end

Base.getindex(output :: SimData, key :: AbstractString) = getindex(output, Symbol(key))
Base.getindex(output :: SimData, inds...)  = getindex(output.data, inds...)
Base.setindex!(output :: SimData, inds...) = setindex!(output.data, inds...)

function Base.show(io :: IO, output :: SimData)
  ids        = tuple(map(string, collect(keys(output.id2ind)))...)
  t, data    = get_data(output)
  d1, d2, d3 = size(data)

  print(io, "SimData{", d1, ",", d2, ",", d3,"}\n")
  print(io, " * species  = ", ids, "\n")
  print(io, " * time     = ", t[end], "\n")
  print(io, " * epochs   = ", d2, "\n")
  print(io, " * trials   = ", d3)
end

@inbounds function update!(
    output    :: SimData,
    Xt        :: Vector,
    t         :: Real,
    epoch     :: Integer,
    trial     :: Integer
)
  t_index, data = get_data(output)
  max_epoch     = length(t_index)

  while (epoch <= max_epoch) && (t >= t_index[epoch])
    update!(data, Xt, epoch, trial)
    epoch = epoch + 1
  end

  return epoch
end

@inline function update!(
    output :: AbstractArray,
    Xt     :: Vector,
    epoch  :: Integer,
    trial  :: Integer
  )
  @inbounds for i in eachindex(Xt)
    output[i, epoch, trial] = Xt[i]
  end
end

"""
```
make_dataframe(x :: SimData)
```

Retrieve simulation data as a `DataFrame`.

Each row represents an observation at a particular `time` and `trial`. Observations for each species are stored in columns.
"""
function get_dataframe(output :: SimData)
    t, data = get_data(output)
    id2ind  = output.id2ind

    d1, d2, d3 = size(data)

    ids = collect(keys(id2ind))

    df = DataFrame()
    df[:time] = repeat(collect(t), outer=[d3])
    for i = 1:d1
        df[ids[i]] = vec(convert(Array, view(data, i, :, :)))
    end
    df[:trial] = repeat(collect(1:d3), inner=[d2])
    return df
end
