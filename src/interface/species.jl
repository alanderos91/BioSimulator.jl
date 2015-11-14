import Base.show

type Species
  id::ASCIIString
  initial::Int
  istracked::Bool
end

function Base.show(io::IO, x::Species)
  if x.istracked
    @printf io "%6d; tracked" x.initial
  else
    @printf io "%6d; untracked" x.initial
  end
end
