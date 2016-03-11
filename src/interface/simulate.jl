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
  stracked::Vector{Int}

  rname::Vector{Symbol}
  rtracked::Vector{Int}

  state::Vector{Int}
  rxns::ReactionVector
  param::Dict{Symbol,Parameter}
end

function Simulation(x::Network)
  state, sname, stracked, inds = _svector(x.species)
  rxns,  rname, rtracked       = _rvector(x.reactions, inds)
  return Simulation(x.id, deepcopy(state), sname, stracked, rname, rtracked, state, rxns, deepcopy(x.parameters))
end

function _svector(species)
  state    = Array(Int,  length(species))
  sname    = Array(Symbol, length(species))
  stracked = Int[]
  inds     = Dict{Symbol,Int}()

  i = 1
  for (key,s) in species
    inds[key] = i
    state[i] = s.population
    sname[i] = s.id
    if s.istracked; push!(stracked, i); end
    i = i + 1
  end
  return state, sname, stracked, inds
end

function _rvector(reactions, inds)
  rxns    = Array(ReactionChannel, length(reactions))
  rname   = Array(Symbol, length(reactions))
  rtracked = Int[]

  j = 1
  for (key,r) in reactions
    rname[j] = r.id
    pre  = Array(Int, length(inds))
    post = Array(Int, length(inds))

    d1 = getfield(r, :reactants)
    d2 = getfield(r, :products)

    for (sname,i) in inds
      pre[i]  = get(d1, sname, 0)
      post[i] = get(d2, sname, 0)
    end
    if r.istracked; push!(rtracked, j); end
    rxns[j] = ReactionChannel(r.id, r.rate, pre, post)
    j = j + 1
  end
  return rxns, rname, rtracked
end

make_observers(::Explicit, sname, stracked, rname, rtracked, spcs, rxns, n, itr) = make_observers(sname, stracked, rname, rtracked, spcs, rxns, 0)

make_observers(::Uniform, sname, stracked, rname, rtracked, spcs, rxns, n, itr) = make_observers(sname, stracked, rname, rtracked, spcs, rxns, n*itr)

make_observers(::Histogram, sname, stracked, rname, rtracked, spcs, rxns, n, itr) = make_observers(sname, stracked, rname, rtracked, spcs, rxns, n)

function make_observers(sname, stracked, rname, rtracked, spcs, rxns, n)
    overseer = Overseer(TimeObserver(:Time, n))

    for i in eachindex(stracked)
        push!(overseer.s_observers, SpeciesObserver(sname[i], spcs, stracked[i], n))
    end

    for i in eachindex(rtracked)
        push!(overseer.r_observers, PropensityObserver(rname[i], rxns[rtracked[i]], n))
    end

    return overseer
end

"""
```
simulate(m::Newtork; with=:ssa, tf=1.0, output=Uniform(), dt=1.0, itr=1, kwargs...)
```

### Arguments
- `m`: The network object

### Optional Arguments
- `with`: The algorithm used in the simulation.
- `tf`: The stopping time for the simulation.
- `output`: The type of output returned by this routine. See [Explicit, Uniform, Histogram].
- `dt`: The step size between updates, if used.
- `itr`: The number of realizations.
- `kwargs`: Optional keyword arguments specific to a particular algorithm. Consult docs for an algorithm for more details.
"""
function simulate(model::Network; with::Symbol=:ssa, tf=1.0, output=Uniform(), dt=1.0, itr=1, kwargs...)
  args = Dict{Symbol,Any}(kwargs)

  if with == :ssa
    algorithm = SSA(args)
  elseif with == :odm
    algorithm = ODM(args)
  elseif with == :frm
    algorithm = FRM(args)
  elseif with == :nrm
    algorithm = NRM(args)
  elseif with == :sal
    algorithm = SAL(args)
  else
    error("$with is an unrecognized algorithm.")
  end
  _run(Simulation(model), algorithm, output, dt, tf, itr)
end

function _run(model::Simulation, alg::Algorithm, output::OutputType, dt, tf, itr)
  # Unpack model
  initial  = model.initial
  spcs     = model.state
  rxns     = model.rxns
  params   = model.param
  sname    = model.sname
  stracked = model.stracked
  rname    = model.rname
  rtracked = model.rtracked

  n  = round(Int, tf / dt) + 1

  overseer = make_observers(output, sname, stracked, rname, rtracked, spcs, rxns, n, itr)
  u  = init_updater(output, overseer, dt, n, itr) # Initialize update manager

  init(alg, rxns, spcs, initial, params)
  for i = 1:itr
    t = 0.0
    copy!(spcs, initial) # Reset copy numbers to initial values
    reset(alg, rxns, spcs, params) # Reset algorithm variables
    while t < tf
      update!(output, u, t)   # Record current state
      τ = step(alg, rxns, spcs, params, t, tf) # Carry out one step of the algorithm
      t = t + τ
    end
    final_update!(output, u, t) # Record final state
  end
  compute_statistics(alg)

  sdata, pdata = compile_data(overseer)
  mdata = compile_metadata(alg, tf, n, itr)
  return SimulationOutput(sdata, pdata, mdata)
end
