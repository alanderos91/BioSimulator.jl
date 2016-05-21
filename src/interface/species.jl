import Base.show

"""
```
Species(id::Symbol, value::Int=0)
```

Construct a Species object named `id` with initial copy number `value`.

### Arguments
- `id`: The name of the `Species`.

### Optional Arguments
- `value`: Initial copy number. Defaults to `0`.
"""
type Species
  id::Symbol
  population::Int

  function Species(id::Symbol, value::Int=0)
    if value < 0
      error("Species population must be nonnegative.")
    end
    new(id, value)
  end
end

function Base.show(io::IO, x::Species)
  print(io, x.population)
end
