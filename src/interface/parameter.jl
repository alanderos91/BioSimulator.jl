import Base.show, Base.convert, Base.promote_rule
import Base.+, Base.-, Base.*, Base./, Base.^, Base.<, Base.>, Base.<=, Base.>=, Base.==
import Base.exp, Base.log, Base.abs, Base.-

"""
```
Parameter(id, value; label)
```

Construct a `Parameter` with identifier `id` that behaves like `Float64(value)`. An optional `label` may be provided as a an annotation.

### Arguments
- `id`: An identifier.
- `value`: A numeric value.
- `label`: A human-readable annotation.
"""
type Parameter
  id    :: Symbol
  value :: Float64
  label :: UTF8String
end

id(p::Parameter) = p.value
value(p::Parameter) = p.value
label(p::Parameter) = p.label

typealias Parameters Dict{Symbol, Parameter}

function Parameter(id::Symbol, value; label="")
  return Parameter(id, Float64(value), label)
end

function Base.show(io::IO, p::Parameter)
  str = (label(p) == "") ? label(p) : string("; ", label(p))
  println(io, id(p), "=", value(p), str)
end

Base.convert{T<:Number}(::Type{T}, p::Parameter) = convert(T, p.value)

Base.promote_rule{T<:Number}(::Type{Parameter}, ::Type{T}) = Float64

for op in [:(Base.(:+)),
           :(Base.(:-)),
           :(Base.(:*)),
           :(Base.(:/)),
           :(Base.(:^)),
           :(Base.max),
           :(Base.min),
           :(Base.(:<)),
           :(Base.(:>)),
           :(Base.(:<=)),
           :(Base.(:>=)),
           :(Base.(:(==)))]

    @eval ($op)(p::Parameter, q::Parameter) = ($op)(p.value, q.value)
    @eval ($op)(p::Parameter, x::Integer)   = ($op)(p.value, x)
    @eval ($op)(p::Parameter, x::Number)    = ($op)(p.value, x)
    @eval ($op)(x::Number, p::Parameter)    = ($op)(x, p.value)
end

for f in [:(Base.exp),
          :(Base.log),
          :(Base.abs),
          :(Base.(:-))]

    @eval ($f)(p::Parameter) = ($f)(p.value)
end
