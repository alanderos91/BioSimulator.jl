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

function distance(
  x :: Site{2,T,Hexagonal2D{T}},
  y :: Site{2,T,Hexagonal2D{T}}) where T
  d1 = x[1] - y[1]
  d2 = x[2] - y[2]

  return sign(d1) == sign(d2) ? abs(d1 + d2) : max(abs(d1), abs(d2))
end

function build_local_neighborhoods!(sites::Vector{Site{2,T,Hexagonal2D{T}}}) where T
  n = length(sites)

  for i in 1:n
    for j in 1:i-1 # j < i
      x = sites[i]
      y = sites[j]
      d = distance(x, y)

      if d == 1
        add_neighbor!(x, y)
        add_neighbor!(y, x)
      end
    end
  end
  return sites
end

### TODO: 3D