immutable SimConfig
  alg_name  :: UTF8String
  sim_time  :: Float64
  wall_time :: Float64
end

immutable SimData
  Xt      :: Array{Int, 3}
  species :: Vector{Symbol}
  id2ind  :: Dict{Symbol,Int}
  #meta    :: SimConfig
end

  function SimData(species, samples, trials, list, dict)
    new(Array(Int, species, samples, trials), list, dict)
  end
end

function Base.show(io::IO, data::SimData)
  m, n, p = size(data.Xt)

  print("SimData[$m species, $n epochs, $p trials]")
end

function getindex(x::SimData, key::Symbol)
  ix = x.id2ind[key]
  x.Xt[ix]
end

function getindex(x::Simdata)
