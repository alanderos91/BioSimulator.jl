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

"""
```
simulate(m::Newtork; with=:ssa, tf=1.0, output=Uniform(), dt=1.0, itr=1, kwargs...)
```

### Arguments
- `m`: The network object

### Optional Arguments
- `with`: The algorithm used in the simulation.
- `tf`: The stopping time for the simulation.
- `output`: The type of output returned by this routine. See [Explicit, Uniform, Mean, Histogram].
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

init_updater(::Explicit,  dt, tracked, n, itr) = Updater(dt, tracked, 0)
init_updater(::Uniform,   dt, tracked, n, itr) = Updater(dt, tracked, n)
init_updater(::Mean,      dt, tracked, n, itr) = Updater(dt, tracked, n)
init_updater(::Histogram, dt, tracked, n, itr) = Updater(dt, tracked, itr)

init_df(::Explicit,  sname, tracked, n, itr) = init_df(sname, tracked, 0)
init_df(::Uniform,   sname, tracked, n, itr) = init_df(sname, tracked, n * itr)
init_df(::Mean,      sname, tracked, n, itr) = init_df(sname, tracked, n)
init_df(::Histogram, sname, tracked, n, itr) = init_df(sname, tracked, itr)

_prepare_df(::OutputType, df, n, itr) = df

function _prepare_df(::Mean, df, n, itr)
  for key in names(df)
    df[key] = df[key] / itr
  end
  return df
end

function _run(model::Simulation, alg::Algorithm, output::OutputType, dt, tf, itr)
  # Unpack model
  initial = model.initial
  sname   = model.sname
  tracked = model.tracked
  spcs    = model.state
  rxns    = model.rxns
  params  = model.param

  n  = round(Int, tf / dt) + 1

  u  = init_updater(output, dt, tracked, n, itr) # Initialize updater
  df = init_df(output, sname, tracked, n, itr)   # Initialize output
  init(alg, rxns, spcs, initial, params)
  for i = 1:itr
    t = 0.0
    copy!(spcs, initial) # Reset copy numbers to initial values
    reset(alg, rxns, spcs, params)    # Reset algorithm variables
    while t < tf
      update!(output, df, t, u, spcs)   # Record current state
      τ = step(alg, rxns, spcs, params, t, tf) # Carry out one step of the algorithm
      t = t + τ
    end
    final_update!(output, df, t, u, spcs) # Record final state
  end
  _prepare_df(output, df, n, itr)
  return df
end
