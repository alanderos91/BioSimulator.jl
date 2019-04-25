struct IPSReactionIR
  pairwise   :: Bool
  input1     :: Int # center
  input2     :: Int # adjacent
  output1    :: Int # center
  output2    :: Int # adjacent
  label      :: Int # reaction index
  p_index    :: Int # parameter index
end

##### helper functions #####
function clean_expression(input_ex)
  ex = Expr(input_ex.head)
  
  for line in input_ex.args
    if line isa Expr && line.head == :tuple && length(line.args) == 2
      # change 0 to :∅ to make it its own "type"
      clean_line = postwalk(x -> x isa Integer && x == 0 ? :∅ : x, line)
      
      # add the new line
      push!(ex.args, clean_line)
    end
  end
  
  # add checks to make sure each formula is well-formed
  
  return ex
end

# catch-all for symbols and constants
function add_tokens!(tokens, sym)
  push!(tokens, sym)
end

# search through an Expr
function add_tokens!(tokens, ex::Expr)
  for token in ex.args
    if isa(token, Symbol)
      push!(tokens, token)
    end
  end
  
  return tokens
end

# build a Dictionary that maps a particle type to an index
function build_species_dict(ex)
  tokens = Symbol[:∅]
  
  for line in ex.args
    formula = line.args[1]
    
    input  = formula.args[1]
    output = formula.args[2]
    
    add_tokens!(tokens, input)
    add_tokens!(tokens, output)
  end
  
  tokens = unique(tokens)
  
  filter!(x -> (x != :+), tokens)
  filter!(x -> (x != :*), tokens)
  
  return OrderedDict{Symbol,Int}(tokens[i] => i for i in eachindex(tokens))
end

# catch-all for non-pairwise events
function get_species_types(sym, dict)
  x = dict[sym]
  y = 0
  
  return x, y
end

# for pairwise events
function get_species_types(ex::Expr, dict)
  x_sym = ex.args[2]
  y_sym = ex.args[3]
  
  x = dict[x_sym]
  y = dict[y_sym]
  
  return x, y
end

# build a Dictionary that maps a parameter symbol to an index
function build_parameters_dict(p)
  parameters = OrderedDict{Symbol,Int}()
  p_count = 0
  for parameter in p
      if !haskey(parameters,parameter)
        p_count += 1
        parameters[parameter] = p_count
      end
  end

  return parameters
end

# build a list of internal IPSReactionIR objects from input
function encode_reaction_struct(ex::Expr, species_dict, params_dict)
  reactions = IPSReactionIR[]
  
  for j in eachindex(ex.args)
    line = ex.args[j]
    
    formula   = line.args[1]
    parameter = line.args[2]
    
    input  = formula.args[1]
    output = formula.args[2]
    
    is_pairwise = isa(input, Expr) # if not, then the argument is a symbol or constant
    
    type1, type2 = get_species_types(input, species_dict)
    type3, type4 = get_species_types(output, species_dict)
    
    reaction = IPSReactionIR(is_pairwise, type1, type2, type3, type4, j, params_dict[parameter])

    push!(reactions, reaction)
  end
  
  return reactions
end

## define a reaction set from user input
macro def_reactions(inputex::Expr, p...)
  __def_reactions(inputex, p)
end

## define a reaction set from a programatically generated list
macro def_reactions(inputex::Symbol, p...)
  escex = esc(inputex)

  quote
    __def_reactions($escex, $p)
  end
end

# implements the body of the @def_reactions macro
function __def_reactions(inputex, p)
  # sweep through the user's code block to clean it up
  # e.g. strip away the annoying line numbers that get folded in
  ex = clean_expression(inputex)
  
  # sweep through the user's code block to pick up
  # all the unique species and map them to an index
  species_dict = build_species_dict(ex)
  params_dict  = build_parameters_dict(p)
  
  # translate the user's model to some intermediate representation
  reactions = encode_reaction_struct(ex, species_dict, params_dict)
  
  # the macro needs to return an expression
  # which builds a IPSReactionIR array
  return :($reactions)
end

## enumerate the full reaction list using a given spatial structure
macro enumerate_with_nclass(r, n, d, p)
  escr = esc(r)
  escn = esc(n)
  escd = esc(d)
  escp = esc(p)
  quote
    __reactions_nclass($escr, $escn, $escd, $escp)
  end
end

# implements the body of @enumerate_with_nclass
function __reactions_nclass(initial, nbtype, d, params)
  nbtype ∉ NBTYPES && error("unsupported neighborhood structure")
  d > 3 && error("do you really need $(d) dimensions?")

  # determine number of particle types
  L = 0

  for reaction in initial
    L = max(L, reaction.input1, reaction.input2, reaction.output1, reaction.output2)
  end

  # determine maximum number of neighbors
  nbmax = 0
  if nbtype == :vonneumann
    nbmax = 2 * d
  elseif nbtype == :hexagon
    d == 1 && (nbmax = 2)
    d == 2 && (nbmax = 6)
    d == 3 && (nbmax = 9)
  end

  # compositions = collect(multiexponents(L + 1, 2 * d))
  compositions = collect(multiexponents(L, nbmax))
  number_compositions = length(compositions)
  
  reactions = IPSReactionStruct[]
  
  for reaction in initial
    for class in 1:number_compositions
      # class is the integer corresponding to each composition
      composition = compositions[class]

      # check for pairwise reaction
      if reaction.pairwise == true
        # get the number of reactants in the composition
        number_reactants = composition[reaction.input2]

        # ignore this class-reaction pair if there are no suitable
        # reactants in the given configuration
        if number_reactants != 0
          rate = params[reaction.p_index]

          push!(reactions, IPSReactionStruct( 
            reaction.pairwise,
            reaction.input1,                      
            reaction.input2,                     
            reaction.output1,                     
            reaction.output2,
            class,                  # neighborhood class
            number_reactants * rate # local rate
          ) )
        end
      else
        rate = params[reaction.p_index]

        push!(reactions, IPSReactionStruct( 
          reaction.pairwise,
          reaction.input1,                      
          reaction.input2,                     
          reaction.output1,                     
          reaction.output2,                      
          class,                     # neighborhood class
          rate / number_compositions # scaled local rate
        ) )
      end
    end
  end

  return reactions
end

## enumerate the full reaction list using a given spatial structure
macro enumerate_with_sclass(r, n, d, p)
  escr = esc(r)
  escn = esc(n)
  escd = esc(d)
  escp = esc(p)
  quote
    __reactions_sclass($escr, $escn, $escd, $escp)
  end
end

# implements the body of the @enumerate_with_sclass macro
function __reactions_sclass(initial, nbtype, d, params)
  nbtype ∉ NBTYPES && error("unsupported neighborhood structure")
  d > 3 && error("do you really need $(d) dimensions?")

  # determine maximum number of neighbors
  nbmax = 0
  if nbtype == :vonneumann
    nbmax = 2 * d
  elseif nbtype == :hexagon
    d == 1 && (nbmax = 2)
    d == 2 && (nbmax = 6)
    d == 3 && (nbmax = 9)
  end

  pairs     = build_reactant_pairs(initial)
  reac2sidx = build_pair_dictionary(pairs, nbmax)
  reactions = IPSReactionStruct[]
  
  for reaction in initial
    sampleidx = reac2sidx[(reaction.input1, reaction.input2)]
    
    if reaction.pairwise == true
      
      # iterate over the possible number of adjacent reactants
      for number_reactants in 1:nbmax
        rate = params[reaction.p_index]

        push!(reactions, IPSReactionStruct( 
        reaction.pairwise,
        reaction.input1,                      
        reaction.input2,                     
        reaction.output1,                     
        reaction.output2,
        # need to shift index based on dimension
        sampleidx + number_reactants - 1,
        # total rate at which a particle in this class undergoes this reaction
        number_reactants * rate)
        )
      end
    else
      # don't care about neighborhood composition here for the rate or the class, so we explicitly ignore it
      # meaning we don't need to change anything about the reaction vector
      rate = params[reaction.p_index]

      push!(reactions, IPSReactionStruct( 
      reaction.pairwise,                    # all as expected...
      reaction.input1,                      
      reaction.input2,                    
      reaction.output1,                     
      reaction.output2,                     
      sampleidx,                              
      rate                          # total rate at which a particle in this class undergoes this reaction
      ) )
    end
  end
  
  return reactions
end
