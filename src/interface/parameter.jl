import Base.show, Base.convert, Base.promote_rule
import Base.+, Base.-, Base.*, Base./, Base.^, Base.<, Base.>, Base.<=, Base.>=, Base.==
import Base.exp, Base.log, Base.abs, Base.-

type Parameter
  id::Symbol
  value::Float64

  description::AbstractString
end

typealias Parameters Dict{Symbol, Parameter}

function parameter(id::Symbol, value::Float64; description="")
  return Parameter(id, value, description)
end

function Base.show(io::IO, p::Parameter)
  @printf io " * %s= %f   (%s)" p.id p.value p.description
end

Base.convert{T<:Number}(::Type{T}, p::Parameter) = convert(T, p.value)

Base.promote_rule{T<:Number}(::Type{Parameter}, ::Type{T}) = promote_rule(Float64, T)

for op in [:(Base.(:+)),
           :(Base.(:-)),
           :(Base.(:*)),
           :(Base.(:/)),
           :(Base.(:^)),
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
