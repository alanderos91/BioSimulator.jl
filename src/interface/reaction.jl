"""
```
Reaction(name, rate, formula)
```

Define a `Reaction` with a `name`, a stochastic `rate` constant, and a `formula`.

The `formula` should use `Species` from a `Network`. Multiple reactants are separated by a `+` symbol; similarly for products. Reactants and products are separated by a `-->` symbol. Refer to the examples below:

# Examples
- Reaction("birth", 1.0, "X --> X + X")
- Reaction("death", 0.5, "X --> 0") # no products
- Reaction("immigration", 0.5, "0 --> X") # no reactants
"""
struct Reaction
  id   :: String
  rate :: Float64

  reactants :: Dict{Symbol,Int}
  products  :: Dict{Symbol,Int}

  function Reaction(name, rate, formula)
    reactants, products = parse_reaction(formula)

    if isempty(reactants) && isempty(products)
      error("Reaction needs at least one reactant or product.")
    end

    if any(x -> x < 0, values(reactants)) || any(x -> x < 0, values(products))
      error("Coefficients must be positive.")
    end
    return new(string(name), rate, reactants, products)
  end

  Reaction(name, rate, reactants, products) = new(string(name), rate, reactants, products)
end

function Base.show(io::IO, x::Reaction)
  print_participants(io, x.reactants)
  print(io, " --> ")
  print_participants(io, x.products)
end

function print_participants(io, participants)
  n = length(participants)
  if n == 0
    print(io, "∅")
  else
    i = 1
    for (id, coeff) in participants
      if coeff == 1
        print(io, id)
      else
        print(io, coeff, " ", id)
      end

      if i < n
        print(io, " + ")
      end
      i = i + 1
    end
  end
end

parse_reaction(formula::String) = parse_reaction(parse(formula))

function parse_reaction(ex::Expr)
  reactants = Dict{Symbol,Int}()
  products  = Dict{Symbol,Int}()

  if ex.head == :-->
    exr = ex.args[1]
    exp = ex.args[2]

    add_participants!(reactants, exr)
    add_participants!(products,  exp)
  else
    throw("malformed reaction")
  end
  return reactants, products
end

function add_participants!(dict, ex)
  if isa(ex, Symbol)
    #push!(dict, ex)
    val = get(dict, ex, 0)
    dict[ex] = val + 1
  elseif isa(ex, Expr)
    if ex.args[1] == :*
      #push!(dict, ex.args[3], ex.args[2])
      val = get(dict, ex.args[3], 0)
      dict[ex.args[3]] = val + ex.args[2]
    else
      for i in 2:length(ex.args)
        add_participants!(dict, ex.args[i])
      end
    end
  end
end
