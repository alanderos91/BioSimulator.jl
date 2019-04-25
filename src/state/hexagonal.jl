struct Hexagonal2D{T} <: Neighborhood
  neighbors :: Vector{Site{2,T,Hexagonal2D{T}}}
  capacity  :: Int

  Hexagonal2D{T}() where T = new{T}(sizehint!(Site{2,T,Hexagonal2D{T}}[], 6), 6)
end

capacity(x :: Type{Hexagonal2D{T}} ) where {T} = 6

@inline function eachdir(::Type{Hexagonal2D{T}}, x) where T
  return ((x[1] - 1, x[2]),
          (x[1] + 1, x[2]),
          (x[1], x[2] - 1),
          (x[1], x[2] + 1),
          (x[1] + 1, x[2] - 1),
          (x[1] - 1, x[2] + 1))
end

### TODO: 3D