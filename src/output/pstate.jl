immutable PopulationState
  tval::Float64
  pop::Int
end

function Base.show(io::IO, s::PopulationState)
  @printf io "%4.4e  %9.0d\n" s.tval s.pop
end
