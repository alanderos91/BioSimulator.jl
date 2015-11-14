import Base.show

type ReactionDef
  id::ASCIIString
  rate::ASCIIString
  reactants::Dict{ASCIIString,Int}
  products::Dict{ASCIIString,Int}

  function ReactionDef(id, rate)
    r = Dict{ASCIIString,Int}()
    p = Dict{ASCIIString,Int}()
    return new(id, rate, r, p)
  end
end

function Base.show(io::IO, x::ReactionDef)
  print_participants(io, x.reactants)
  print(io, " ---> ")
  print_participants(io, x.products)
end

function print_participants(io, participants)
  n = length(participants)
  if n == 0
    print("∅")
  else
    i = 1;
    for (s,c) in participants
      print(io,c," ",s)
      if i < n; print(io," + "); end
      i = i + 1
    end
  end
end

function add_participant!(r, sname, coeff, typ)
  d = if typ == :reactants
    getfield(r, typ)
  else
    getfield(r, :products)
  end
  setindex!(d, coeff, sname)
end

function rmv_participant!(r, sname, typ)
  d = if typ == :reactants
    getfield(r, typ)
  else
    getfield(r, :products)
  end
  delete!(d, sname)
end
