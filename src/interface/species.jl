import Base.show

"""
```
Species(id, [value=0])
```

Define a `Species` with a `name` and initial population `value`.
"""
type Species
  id         :: String
  population :: Int

  function Species(id, value::Int=0)
    if value < 0
      error("Species population must be nonnegative.")
    end
    new(string(id), value)
  end
end

id(x::Species) = x.id
population(x::Species) = x.population

function Base.show(io::IO, x::Species)
  println(io, x.population)
end
