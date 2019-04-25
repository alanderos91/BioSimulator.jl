struct VonNeumann1D{T} <: Neighborhood
  neighbors :: Vector{Site{1,T,VonNeumann1D{T}}}
  capacity  :: Int

  VonNeumann1D{T}() where T = new(sizehint!(Site{1,T,VonNeumann1D{T}}[], 2), 2)
end

function eachdir(::Type{VonNeumann1D{T}}, x) where T
  return ((x[1] - 1,),
          (x[1] + 1,))
end

struct VonNeumann2D{T} <: Neighborhood
  neighbors :: Vector{Site{2,T,VonNeumann2D{T}}}
  capacity  :: Int

  VonNeumann2D{T}() where T = new(sizehint!(Site{2,T,VonNeumann2D{T}}[], 4), 4)
end

@inline function eachdir(::Type{VonNeumann2D{T}}, x) where T
  return ((x[1] - 1, x[2]),
          (x[1] + 1, x[2]),
          (x[1], x[2] - 1),
          (x[1], x[2] + 1))
end

struct VonNeumann3D{T} <: Neighborhood
  neighbors :: Vector{Site{3,T,VonNeumann3D{T}}}
  capacity  :: Int

  VonNeumann3D{T}() where T = new(sizehint!(Site{3,T,VonNeumann3D{T}}[], 6), 6)
end

@inline function eachdir(::Type{VonNeumann3D{T}}, x) where T
  return ((x[1] - 1, x[2], x[3]),
          (x[1] + 1, x[2], x[3]),
          (x[1], x[2] - 1, x[3]),
          (x[1], x[2] + 1, x[3]),
          (x[1], x[2], x[3] - 1),
          (x[1], x[2], x[3] + 1))
end

capacity(x :: Type{VonNeumann1D{T}} ) where {T} = 2
capacity(x :: Type{VonNeumann2D{T}} ) where {T} = 4
capacity(x :: Type{VonNeumann3D{T}} ) where {T} = 6
