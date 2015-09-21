export Species, reset!

"""
A key-value pair that represents a population.

### Constructors

```julia
Species(id::ASCIIString, pop::Int, istracked::Bool)
```

### Arguments
* `id` : an ASCIIString identifier for the population
* `pop` : an Int representing the population size
* `istracked` : a Bool specifying whether to include in simulation output
"""
type Species
    id::ASCIIString
    pop::Int
    istracked::Bool
end

"""
```reset!(x::Species, y::Species)```

Reset the fields of a Species instance to some initial values.
Note that this does *not* create a new instance.

### Arguments
* `x` : The instance that will hold the copied fields.
* `y` : The instance whose fields will be copied.
"""

function reset!(x::Species, y::Species)
  x.id = y.id
  x.pop = y.pop
  x.istracked = y.istracked
  return x
end

"""
```reset!(dest::Vector{Species}, src::Vector{Species})```

Reset the fields of each Species in a Vector{Species} to some initial values.
Same as invoking reset! on each pairing of the two vectors. Vectors must be the same length.
### Arguments
* `dest` : The Vector{Species} that will contain the copied fields.
* `src` : The Vector{Species} representing the initial values.

Based on Base.copy!
"""
function reset!(dest::Vector{Species}, src::Vector{Species})
  n = length(src)
  n > length(dest) && throw(BoundsError(dest, n))
  @inbounds for i = 1:n
    reset!(dest[i], src[i])
  end
  return dest
end
