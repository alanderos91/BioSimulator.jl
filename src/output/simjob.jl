immutable SimulationJob
  results::Vector{SimulationResult}
end

function SimulationJob()
  return SimulationJob(SimulationResult[])
end

function SimulationJob(n::Integer)
  return SimulationJob(Array{SimulationResult}(n))
end

push!(sj::SimulationJob, r::SimulationResult) = push!(sj.results, r)
getindex(sj::SimulationJob, i::Integer) = getindex(sj.results, i)

function setindex!(sj::SimulationJob, r::SimulationResult, i::Integer)
  setindex!(sj.results, r, i)
end

function isempty(sj::SimulationJob)
  return isempty(sj.results)
end

function length(sj::SimulationJob)
  return length(sj.results)
end
