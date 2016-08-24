

type PropensityVector{T <: AbstractFloat} <: AbstractVector{T}
  cache       :: Vector{T}
  intensity   :: T
  error_bound :: T

  PropensityVector(x::Vector{T}) = new(x, sum(x), zero(T))
end

typealias PVec{T} PropensityVector{T}

##### constructor #####
call{T}(::Type{PVec{T}}, m::Integer) = PVec{T}(zeros(T, m))
call{T}(::Type{PVec}, x::Vector{T})  = PVec{T}(x)

##### PVec interface #####
intensity(x::PVec) = x.intensity
islossy(x::PVec)   = eps(eltype(x)) * x.intensity < x.error_bound# for Float64, eps is 2^-32

##### AbstractArray interface #####

Base.size(x::PVec)    = size(x.cache)
Base.size(x::PVec, d) = size(x.cache, d)

Base.getindex(x::PVec, i::Int) = x.cache[i]

@inline function Base.setindex!(x::PVec, v, i::Int)
  # update error bound first
  x.error_bound = eps(eltype(x)) * (x.intensity + x[i] + v)

  # update the cumulative sum
  x.intensity = x.intensity - x[i] + v

  # update the cache
  x.cache[i] = v
end

Base.linearindexing(x::PVec) = Base.LinearFast()
