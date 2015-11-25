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
    sname::Vector{Symbol}
    tracked::Vector{Int}
    inds::Dict{ASCIIString,Int}

    state::Vector{Int}
    rxns::Vector{Reaction}
    param::Dict{ASCIIString,Float64}
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
  inds = Dict{ASCIIString,Int}()

  i = 1
  for (key,s) in species
    inds[key] = i
    state[i] = s.initial
    sname[i] = symbol(s.id)
    if s.istracked; push!(tracked, i); end
    i = i + 1
  end
  return state, sname, tracked, inds
end

function _rvector(reactions, inds)
  rxns = Array(Reaction, length(reactions))

  j = 1
  for (key,r) in reactions
    pre = Array(Int, length(inds))
    post = Array(Int, length(inds))
    for (sname,i) in inds
      pre[i]  = get(r.reactants, sname, 0)
      post[i] = get(r.products,  sname, 0)
    end
    rxns[j] = Reaction(r.id, r.rate, pre, post)
    j = j + 1
  end
  return rxns
end

function simulate(model::Network; with::Symbol=:ssa, tf=1.0, output=Uniform(), dt=1.0, itr=1, kwargs...)
  args = Dict{Symbol,Any}(kwargs)
  
  if with == :ssa
    algorithm = SSA(itr, tf, dt, args)
  elseif with == :odm
    algorithm = ODM(itr, tf, dt, args)
  elseif with == :frm
    algorithm = FRM(itr, tf, dt, args)
  elseif with == :nrm
    algorithm = NRM(itr, tf, dt, args)
  elseif with == :sal
    algorithm = SAL(itr, tf, dt, args)
  else
    error("$with is an unrecognized algorithm.")
  end

  _run(Simulation(model), algorithm, output)
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

function _run(model::Simulation, alg::Algorithm, output::OutputType)
  # Unpack model
  initial = model.initial
  sname   = model.sname
  tracked = model.tracked
  spcs    = model.state
  rxns    = model.rxns
  params  = model.param

  n  = round(Int, alg.tf / alg.dt) + 1

  u  = init_updater(output, alg.dt, tracked, n, alg.itr) # Initialize updater
  df = init_df(output, sname, tracked, n, alg.itr)   # Initialize output
  init(alg, rxns, spcs, initial, params)
  for i = 1:alg.itr
    copy!(spcs, initial) # Reset copy numbers to initial values
    reset(alg, rxns, spcs, params)    # Reset algorithm variables
    while alg.t < alg.tf
      update!(output, df, alg.t, u, spcs)   # Record current state
      step(alg, rxns, spcs, params) # Carry out one step of the algorithm
    end
    final_update!(output, df, alg.t, u, spcs) # Record final state
  end
  _prepare_df(output, df, n, alg.itr)
  return df
end
