import Base.show

type Species
  id::Symbol
  population::Int
  istracked::Bool

  function Species(id, value, istracked)
    if value < 0
      error("Species population must be nonnegative.")
    end
    new(id, value, istracked)
  end
end

"""
```
Species(id, value::Int=0, istracked::Bool=true)
```

Construct a Species object named `id` with initial copy number `value`.

### Arguments
- `id`: The name of the Species.

### Optional Arguments
- `value`: Initial copy number.
- `istracked`: Indicates whether a simulation should record this `Species`'s population.
  - `true`: Track this `Species`.
  - `false`: Do not track this `Species`.
"""
function Species(id, value::Int=0; istracked::Bool=true)
  Species(id, value, istracked)
end

function Base.show(io::IO, x::Species)
  if x.istracked
    @printf io "%6d; tracked"   x.population
  else
    @printf io "%6d; untracked" x.population
  end
end
