import Base.show

type Reaction
  id::Symbol
  rate::Symbol
  reactants::Dict{Symbol,Int}
  products::Dict{Symbol,Int}

  function Reaction(id, rate; r=(), p=())
    reactants = Dict{Symbol,Int}(r)
    products  = Dict{Symbol,Int}(p)

    if isempty(reactants) && isempty(products)
      error("Reaction needs at least one reactant or product.")
    end

    if any(x -> x < 0, values(reactants)) || any(x -> x < 0, values(products))
      error("Coefficients must be positive.")
    end
    return new(id, rate, Dict(reactants), Dict(products))
  end
end

function Base.show(io::IO, x::Reaction)
  print_participants(io, x.reactants)
  @printf io " ---> "
  print_participants(io, x.products)
end

function print_participants(io, participants)
  n = length(participants)
  if n == 0
    @printf io "∅"
  else
    i = 1
    for (id, coeff) in participants
      if coeff == 1
        @printf io "%s" id
      else
        @printf io "%d %s" coeff id
      end

      if i < n
        @printf io " + "
      end
      i = i + 1
    end
  end
end
