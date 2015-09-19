export Species, reset!

type Species
    id::ASCIIString
    pop::Int
    istracked::Bool
end

# Similar to Base.copy!; used to reset populations to initial values.
function reset!(dest::Vector{Species}, src::Vector{Species})
  n = length(src)
  n > length(dest) && throw(BoundsError(dest, n))
  @inbounds for i = 1:n
    dest[i].id = src[i].id
    dest[i].pop = src[i].pop
    dest[i].istracked = src[i].istracked
  end
  return dest
end
