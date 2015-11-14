
  type Reaction
    id::ASCIIString
    rate::ASCIIString
    propensity::Float64
    pre::Vector{Int}
    post::Vector{Int}

    function Reaction(id::ASCIIString, rate::ASCIIString, pre::Vector{Int}, post::Vector{Int})
      if any(pre .< 0) || any(post .< 0)
        error("Stoichiometric coefficients must be positive.")
      end
      return new(id, rate, 0.0, pre, post)
    end
  end

  immutable Simulation
      id::ASCIIString
      initial::Vector{Int}
      tracked::Dict{ASCIIString,Int}
      inds::Dict{ASCIIString,Int}

      state::Vector{Int}
      rxns::Vector{Reaction}
      param::Dict{ASCIIString, Float64}
  end

  function Simulation(x::Network)
    state, tracked, inds = _svector(x.species)
    rxns = _rvector(x.reactions, inds)
    return Simulation(x.id, deepcopy(state), tracked, inds, state, rxns, deepcopy(x.parameters))
  end

  function _svector(species)
    state = Array(Int,  length(species))
    tracked = Dict{ASCIIString,Int}()
    inds = Dict{ASCIIString,Int}()

    i = 1
    for (key,s) in species
      inds[key] = i
      state[i] = s.initial
      if s.istracked; tracked[key] = i; end
      i = i + 1
    end
    return state, tracked, inds
  end

  function _rvector(reactions, inds)
    rxns = Array(Reaction, length(reactions))

    j = 1
    for (key,r) in reactions
      pre = Array(Int, length(inds))
      post = Array(Int, length(inds))
      for (sname,i) in inds
        if haskey(r.reactants, sname)
          pre[i] = r.reactants[sname]
        else
          pre[i] = 0
        end

        if haskey(r.products, sname)
          post[i] = r.products[sname]
        else
          post[i] = 0
        end
      end
      rxns[j] = Reaction(r.id, r.rate, pre, post)
      j = j + 1
    end
    return rxns
  end

  function simulate(model::Simulation, t_final::Float64, with::Symbol=:ssa; o::OutputType=Uniform(), dt::Float64=1.0, itr::Int=1, args...)

  	d = Dict{Symbol, Float64}()
  	for (key, value) in args
  		d[key] = value
  	end

  	if with == :ssa
  		ssa(model, t_final, o, dt, itr)
  	elseif with == :odm
  		c       = haskey(d, :c)       ? d[:c]       : Tight()
  		steps   = haskey(d, :steps)   ? d[:steps]   : 100
  		samples = haskey(d, :samples) ? d[:samples] : 1

  		odm(model, t_final, o, dt, itr, c, steps, samples)
  	elseif with == :frm
  		frm(model, t_final, o, dt, itr)
  	elseif with == :nrm
  		nrm(model, t_final, o, dt, itr)
  	elseif with == :sal
  		tol   = haskey(d, :tol)   ? d[:tol]   : 0.125
  		thrsh = haskey(d, :thrsh) ? d[:thrsh] : 100.0
  		ctrct = haskey(d, :ctrct) ? d[:ctrct] : 0.75

  		sal(model, t_final, o, dt, itr, tol, thrsh, ctrct)
  	end
  end

  export Simulation, simulate
