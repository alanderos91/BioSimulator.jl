import Base.show

"""
```
Reaction(id, rate, r=(), p=())
```

### Arguments
- `id`: An identifier for this `Reaction`.
- `rate`: An identifier for this `Reaction`'s rate constant.

### Optional Arguments
- `r`: A list of reactants and their coefficients, specified using a `Tuple` of `Pair`s e.g. `(:X1 => 1, :X2 => 1)`.
- `p`: Like `r`, but a list for reaction products.
"""
type Reaction
  id::Symbol
  rate::Symbol
  reactants::Dict{Symbol,Int}
  products::Dict{Symbol,Int}

  function Reaction(id, rate, expression::Expr)
    reactants, products = parse_reaction(expression)

    if isempty(reactants) && isempty(products)
      error("Reaction needs at least one reactant or product.")
    end

    if any(x -> x < 0, values(reactants)) || any(x -> x < 0, values(products))
      error("Coefficients must be positive.")
    end
    return new(id, rate, reactants, products)
  end
end

function Base.show(io::IO, x::Reaction)
  print_participants(io, x.reactants)
  @printf io " --> "
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
