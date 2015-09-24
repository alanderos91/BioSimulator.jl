immutable SimulationJob
  results::Vector{SimulationResult}
end

function SimulationJob()
  return SimulationJob(SimulationResult[])
end

push!(sj::SimulationJob, r::SimulationResult) = push!(sj.results, r)
getindex(sj::SimulationJob, i::Integer) = getindex(sj.states, i)

function setindex!(sj::SimulationJob, r::SimulationResult, i::Integer)
  setindex!(sj.states, r, i)
end

function isempty(sj::SimulationJob)
  return isempty(sj.states)
end

function length(sj::SimulationJob)
  return length(sj.states)
end
