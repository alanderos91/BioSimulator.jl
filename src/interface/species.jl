import Base.show

"""
```
Species(id, value::Int=0)
```

Construct a Species object named `id` with initial copy number `value`.

### Arguments
- `id`: The name of the `Species`.

### Optional Arguments
- `value`: Initial copy number. Defaults to `0`.
"""
type Species
  id         :: UTF8String
  population :: Int

  function Species(id, value::Int=0)
    if value < 0
      error("Species population must be nonnegative.")
    end
    new(UTF8String(string(id)), value)
  end
end

id(x::Species) = x.id
value(x::Species) = x.population

function Base.show(io::IO, x::Species)
  println(io, x.population)
end
