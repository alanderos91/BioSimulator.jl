struct SimulationSummary{T}
  model :: Network

  algorithm_name   :: String
  algorithm_params :: Dict{Symbol,Any}
  algorithm_stats  :: Dict{Symbol,Int}

  simulation_data :: T

  id2index :: Dict{Symbol,Int}
end

function SimulationSummary(model, algname, algorithm, time, epochs, trials, id2index, output; kwargs...)
  algorithm_name = "$(algname)"
  algorithm_params = Dict{Symbol,Any}(kwargs)
  algorithm_stats = algorithm.stats

  algorithm_params[:time]   = time
  algorithm_params[:epochs] = epochs
  algorithm_params[:trials] = trials

  simulation_data = output

  return SimulationSummary(
    deepcopy(model),
    algorithm_name,
    algorithm_params,
    algorithm_stats,
    simulation_data,
    id2index
  )
end

function Base.show(io :: IO, result :: SimulationSummary)
  show(io, result.model)

  println(io, "")
  print(io, "[ $(result.algorithm_name) ]\n")
  
  print(io, "  parameters:\n")

  for (key, val) in result.algorithm_params
    print(io, "    ", key, " => ", val, "\n")
  end

  print(io, "  statistics:\n")

  for (key, val) in result.algorithm_stats
    print(io, "    ", key, " => ", val, "\n")
  end
end

get_regular_path(xw :: RegularPath{T1,T2}, t, epochs) where {T1,T2} = xw

function get_regular_path(xw::SamplePath{T1,T2}, t, epochs) where {T1,T2}
  max_epoch = epochs + 1
  tdata = collect(linspace(0.0, t, max_epoch))
  xdata = xw.xdata
  epoch = 1

  while (epoch <= max_epoch) && (t >= tdata[epoch])
    for i in eachindex(x)
      @inbounds xdata[i, epoch] = x[i]
    end
    epoch = epoch + 1
  end
  
  return SamplePath(tdata, xdata)
end

@recipe function f(result :: SimulationSummary{T};
    plot_type = :trajectory,
    trial = nothing,
    species = nothing,
    time = nothing, # unused for now
    epochs = nothing
  ) where T

  if species == nothing
    species = collect(keys(result.id2index))
    labels = map(String, keys(result.id2index))
  else
    labels = map(s -> String(s), species)
  end

  label --> reshape(labels, 1, length(labels))

  if plot_type == :trajectory
    if trial == nothing
      trial = 1
    end

    title --> "trial = $(trial)"

    _trajectory(result, trial, species, time)
  elseif plot_type == :meantrajectory
    if epochs == nothing 
      epochs = 1_000
    end
    time = result.algorithm_params[:time]

    title --> "mean trajectory, n = $(length(result.simulation_data))"

    _meantrajectory(result, species, time, epochs)
  elseif plot_type == :histogram
    # if time == nothing
    #   time = result.algorithm_params[:time]
    # end
    time = result.algorithm_params[:time]

    title --> "histogram at t = $(time)"

    _histogram(result, species, time)
  else
    error("unrecognized plot type $(plot_type)")
  end
end

function _trajectory(result, trial, species, time)
  xw = result.simulation_data[trial]
  idxs = select_species(result, species)
  xw = extract_paths(xw, idxs)

  return xw
end

function _meantrajectory(result, species, time, epochs)
  idxs = select_species(result, species)
  xw = extract_paths(result.simulation_data, idxs)
  yw = [ get_regular_path(xw[i], time, epochs) for i in 1:length(xw) ]

  return AveragePath(yw)
end

function _histogram(result, species, time)
  idxs = select_species(result, species)
  xw = extract_paths(result.simulation_data, idxs)

  return BioSimHistogram(xw)
end

function select_species(result :: SimulationSummary, species)
  id2index = result.id2index
  species_keys = [Symbol(s) for s in species]
  idxs = [id2index[s] for s in species_keys]
end

function extract_paths(xw, idxs)
  return [extract_paths(xw[i], idxs) for i in 1:length(xw)]
end

function extract_paths(xw :: RegularPath, idxs)
  RegularPath(xw.tdata, xw.xdata[idxs, :])
end

function extract_paths(xw :: SamplePath, idxs)
  SamplePath(xw.tdata, map(x -> x[idxs], xw.xdata))
end

function DataFrame(result :: SimulationSummary)
  trials  = result.algorithm_params[:trials]
  sample  = result.simulation_data
  nsteps  = [length(sample[i].tdata) for i in 1:trials]
  species = collect(keys(result.id2index))
  idxs    = collect(values(result.id2index))

  # create columns
  tcol = vcat((sample[i].tdata for i in 1:trials)...)
  xcol = reshape_xdata(sample, trials)
  ncol = [i for i in 1:trials for j in 1:nsteps[i]]

  # build the dataframe
  df = DataFrame()

  # add trial label
  df[:trial] = ncol

  # add time points
  df[:time]  = tcol

  # add state data
  for i in eachindex(species)
    df[species[i]] = xcol[:, idxs[i]]
  end

  return df
end

function reshape_xdata(sample :: RegularEnsemble{T1,T2}, trials) where {T1,T2}
  transpose(hcat((sample[i].xdata for i in 1:trials)...))
end

function reshape_xdata(sample :: Ensemble{T1,T2}, trials) where {T1,T2}
  return transpose(hcat((hcat(sample[i].xdata...) for i in 1:trials)...))
end