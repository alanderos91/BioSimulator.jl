struct VonNeumann end

capacity(::VonNeumann, d::Int) = 2 * d

@inline function eachdir(::VonNeumann, x::Site{1})
  return ((x[1] - 1,),
          (x[1] + 1,))
end

@inline function eachdir(::VonNeumann, x::Site{2})
  return ((x[1] - 1, x[2]),
          (x[1] + 1, x[2]),
          (x[1], x[2] - 1),
          (x[1], x[2] + 1))
end

@inline function eachdir(::VonNeumann, x::Site{3})
  return ((x[1] - 1, x[2], x[3]),
          (x[1] + 1, x[2], x[3]),
          (x[1], x[2] - 1, x[3]),
          (x[1], x[2] + 1, x[3]),
          (x[1], x[2], x[3] - 1),
          (x[1], x[2], x[3] + 1))
end

##### distance - this could be done with @generated instead
@inline function distance(::VonNeumann, x::Site{1}, y::Site{1})
  return abs(x[1] - y[1])
end

@inline function distance(::VonNeumann, x::Site{2}, y::Site{2})
  return abs(x[1] - y[1]) + abs(x[2] - y[2])
end

@inline function distance(::VonNeumann, x::Site{3}, y::Site{3})
  return abs(x[1] - y[1]) + abs(x[2] - y[2]) + abs(x[3] - y[3])
end
