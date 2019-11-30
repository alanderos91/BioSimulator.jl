"""
```
Species(id, [value=0])
```

Define a `Species` with a `name` and initial population `value`.
"""
struct Species
  identifier :: Symbol
  population :: Int

  function Species(id, value::Int=0)
    if value < 0
      error("Species population must be nonnegative.")
    end
    new(Symbol(id), value)
  end
end

function Base.show(io::IO, x::Species)
  println(io, x.population)
end
