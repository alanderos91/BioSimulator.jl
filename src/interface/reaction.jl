import Base.show

"""
```
Reaction(id, k::Float64, formula::Expr)
```

Construct a `Reaction` with identifier `id` and stochastic rate constant `k`, represented by a `formula` (e.g. X + Y --> XY).

### Arguments
- `id`: An identifier for this `Reaction`.
- `rate`: A real-valued rate constant.
- `formula`: A chemical formula representation of the `Reaction`. For example, `:(X --> X + X)` defines a birth reaction with reactant `X`, coefficient 1, and product `X`, coefficient 2; :(X --> 2*X) is also valid syntax. The `-->` symbol must be used to separate reactants and products.
"""
type Reaction
  id   :: UTF8String
  rate :: Float64

  reactants :: Dict{Symbol,Int}
  products  :: Dict{Symbol,Int}

  function Reaction(id, rate, expression::Expr)
    reactants, products = parse_reaction(expression)

    if isempty(reactants) && isempty(products)
      error("Reaction needs at least one reactant or product.")
    end

    if any(x -> x < 0, values(reactants)) || any(x -> x < 0, values(products))
      error("Coefficients must be positive.")
    end
    return new(UTF8String(string(id)), rate, reactants, products)
  end
end

function Reaction(id, rate, str::AbstractString) Reaction(id, rate, parse(str))
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

parse_reaction(str::AbstractString) = parse_reaction(parse(str))

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
