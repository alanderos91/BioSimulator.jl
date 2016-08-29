const TRUNCATION_CONST = 2.0e-32

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
islossy(x::PVec)   = x.error_bound > TRUNCATION_CONST * x.intensity

@fastmath function update_errorbound!{T}(x::PVec{T}, xi::T, i::Integer)
  x.error_bound = x.error_bound + eps(T) * (x.intensity + x[i] + xi)
end

@fastmath function update_intensity!{T}(x::PVec{T}, xi::T, i::Integer)
  x.intensity = x.intensity - x[i] + xi
end

##### AbstractArray interface #####

Base.size(x::PVec)             = size(x.cache)
Base.size(x::PVec, d::Integer) = size(x.cache, d)

Base.getindex(x::PVec, i::Integer) = x.cache[i]

Base.setindex!(x::PVec, xi, i::Integer) = (x.cache[i] = xi)

Base.linearindexing(x::PVec) = Base.LinearFast()
