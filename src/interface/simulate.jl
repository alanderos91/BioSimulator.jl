
  type ReactionChannel
    id::Symbol
    rate::Symbol
    propensity::Float64
    pre::Vector{Int}
    post::Vector{Int}

    function ReactionChannel(id, rate, pre, post)
      if any(pre .< 0) || any(post .< 0)
        error("Stoichiometric coefficients must be positive.")
      end
      return new(id, rate, 0.0, pre, post)
    end
  end

  typealias ReactionVector Vector{ReactionChannel}

  immutable Simulation
      id::ASCIIString
      initial::Vector{Int}
      sname::Vector{Symbol}
      tracked::Vector{Int}
      inds::Dict{Symbol,Int}

      state::Vector{Int}
      rxns::ReactionVector
      param::Dict{Symbol,Parameter}
  end

  function Simulation(x::Network)
    state, sname, tracked, inds = _svector(x.species)
    rxns = _rvector(x.reactions, inds)
    return Simulation(x.id, deepcopy(state), sname, tracked, inds, state, rxns, deepcopy(x.parameters))
  end

  function _svector(species)
    state = Array(Int,  length(species))
    sname = Array(Symbol, length(species))
    tracked = Int[]
    inds = Dict{Symbol,Int}()

    i = 1
    for (key,s) in species
      inds[key] = i
      state[i] = s.population
      sname[i] = s.id
      if s.istracked; push!(tracked, i); end
      i = i + 1
    end
    return state, sname, tracked, inds
  end

  function _rvector(reactions, inds)
    rxns = Array(ReactionChannel, length(reactions))

    j = 1
    for (key,r) in reactions
      pre  = Array(Int, length(inds))
      post = Array(Int, length(inds))

      d1 = getfield(r, :reactants)
      d2 = getfield(r, :products)

      for (sname,i) in inds
        pre[i]  = get(d1, sname, 0)
        post[i] = get(d2, sname, 0)
      end

      rxns[j] = ReactionChannel(r.id, r.rate, pre, post)
      j = j + 1
    end
    return rxns
  end
  # TODO: Move argument parsing within algorithms themselves.
  function simulate(model::Simulation, t_final::Float64, with::Symbol=:ssa; o::OutputType=Uniform(), dt::Float64=1.0, itr::Int=1, kwargs...)
    args = Dict{Symbol,Any}(kwargs)

  	if with == :ssa
  		ssa(model, t_final, o, dt, itr)
  	elseif with == :odm
  		c       = haskey(args, :c)       ? args[:c]       : Tight()
  		steps   = haskey(args, :steps)   ? args[:steps]   : 100
  		samples = haskey(args, :samples) ? args[:samples] : 1

  		odm(model, t_final, o, dt, itr, c, steps, samples)
  	elseif with == :frm
  		frm(model, t_final, o, dt, itr)
  	elseif with == :nrm
  		nrm(model, t_final, o, dt, itr)
  	elseif with == :sal
  		tol   = haskey(args, :tol)   ? args[:tol]   : 0.125
  		thrsh = haskey(args, :thrsh) ? args[:thrsh] : 100.0
  		ctrct = haskey(args, :ctrct) ? args[:ctrct] : 0.75

  		sal(model, t_final, o, dt, itr, tol, thrsh, ctrct)
  	end
  end

  export Simulation, simulate
