struct Hexagonal end

capacity(::Hexagonal, d::Int) = 3 * d

@inline function eachdir(::Hexagonal, x::Site{2})
  return ((x[1] - 1, x[2]),
          (x[1] + 1, x[2]),
          (x[1], x[2] - 1),
          (x[1], x[2] + 1),
          (x[1] + 1, x[2] - 1),
          (x[1] - 1, x[2] + 1))
end

function distance(::Hexagonal, x::Site{2}, y::Site{2})
  d1 = x[1] - y[1]
  d2 = x[2] - y[2]

  return sign(d1) == sign(d2) ? abs(d1 + d2) : max(abs(d1), abs(d2))
end
