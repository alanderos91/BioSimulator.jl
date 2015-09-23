export PopulationTrace, SimulationResult, plot, regularize, DataFrame

immutable PopulationState
  name::ASCIIString
  time::Float64
  value::Int
end

function Base.show(io::IO, s::PopulationState)
  @printf io " %8s   %4.4e  %9.0d\n" string(s.name) s.time s.value
end



immutable PopulationTrace
  states::Vector{PopulationState}
end

function PopulationTrace()
  return PopulationTrace(PopulationTrace[])
end

function PopulationTrace(x::Species)
  state = PopulationState(x.id, 0.0, x.pop)
  return PopulationTrace([state])
end

Base.push!(t::PopulationTrace, s::PopulationState) = push!(t.states, s)
Base.getindex(t::PopulationTrace, i::Integer) = getindex(t.states, i)

function Base.setindex!(t::PopulationTrace, s::PopulationState, i::Integer)
  setindex!(t.states, s, i)
end

function Base.isempty(t::PopulationTrace)
  return isempty(t.states)
end

function Base.length(t::PopulationTrace)
  return length(t.states)
end

function species(t::PopulationTrace)
  return t.states[1].name
end

function toarrays(tr::PopulationTrace)
  tval = Float64[]
  pval = Int[]

  for state in tr.states
    push!(tval, state.time)
    push!(pval, state.value)
  end
  return tval, pval
end

function Base.show(io::IO, t::PopulationTrace)
  @printf io "   Name        Time       Value  \n"
  @printf io "----------  ----------  ---------\n"
  for state in t.states
    show(io, state)
  end
end

function init_traces(spcs::Vector{Species})
  traces = Dict{ASCIIString, PopulationTrace}()
  for s in spcs
    if s.istracked
      traces[s.id] = PopulationTrace(s)
    end
  end
  return traces
end

function update_traces!(traces::Dict{ASCIIString, PopulationTrace}, t::Float64, spcs::Vector{Species}, store_trace::Bool)
  for s in spcs
    s.istracked ? update!(traces[s.id], t, s, store_trace) : Nothing
  end
end

function update!(tr::PopulationTrace, t::Float64, s::Species, store_trace::Bool)
  state = PopulationState(s.id, t, s.pop)

  if store_trace
    push!(tr, state)
  end
end

function regularize(tr::PopulationTrace, stepsize::Float64, t_final::Float64)
  ntr = PopulationTrace()

  t = 0.0
  s = tr[1]
  count = 0
  maxcount = round(Int, t_final / stepsize) + 1
  i = 1

  while count < maxcount
    if s.time >= t
      push!(ntr, PopulationState(s.name, t, s.value))
      t = t + stepsize
      count = count + 1
    else
      i = i+1
      s = tr[i]
    end
  end

  return ntr
end

import Gadfly.ElementOrFunctionOrLayers

function plot(tr::PopulationTrace, elements::ElementOrFunctionOrLayers...)
  if !isempty(tr)
    name = species(tr)
    tval, pval = toarrays(tr)
    return Gadfly.plot(x=tval, y=pval, elements...)
  else
    # TODO
  end
end

function DataFrame(tr::PopulationTrace; itr::Int=1)
  if !isempty(tr)
    s = species(tr)
    tval, pval = toarrays(tr)
    return DataFrames.DataFrame(Time=tval, Population=pval, Species=s, Realization=itr)
  else
    # TODO
  end
end

function DataFrame(traces::Dict{ASCIIString,PopulationTrace}; itr::Int=1)
  if !isempty(traces)
    df = DataFrames.DataFrame()
    for tr in values(traces)
      df = vcat(df, DataFrame(tr, itr=itr))
    end
    return df
  else
    # TODO
  end
end

immutable SimulationResults
  model::ASCIIString
  results::Vector{Dict{ASCIIString,PopulationTrace}}
  metadata::Dict{Any,Any}
end

function Base.show(io::IO, sr::SimulationResults)
  @printf io "[ %s ]\n" sr.model
  for (key,val) in sr.metadata
    @printf io " * %s: %s\n" key val
  end
end

function regularize(trs::Dict{ASCIIString,PopulationTrace}, stepsize::Float64, t_final::Float64)
  ntrs = Dict{ASCIIString, PopulationTrace}()
  for (key, trace) in trs
    ntrs[key] = regularize(trace, stepsize, t_final)
  end
  return ntrs
end

function regularize(rslt::Vector{Dict{ASCIIString,PopulationTrace}}, stepsize::Float64, t_final::Float64)
  results = Dict{ASCIIString,PopulationTrace}[]
  for r in rslt
    d = regularize(r, stepsize, t_final)
    push!(results, d)
  end
  return results
end

function regularize(sr::SimulationResults, stepsize::Float64, t_final::Float64)
  alg = sr.model
  md = sr.metadata
  results = regularize(sr.results, stepsize, t_final)
  return SimulationResults(alg, results, md)
end

function DataFrame(sr::SimulationResults)
  df = DataFrames.DataFrame()
  for i in eachindex(sr.results)
    df = vcat(df, DataFrame(sr.results[i], itr=i))
  end
  return df
end
