immutable SimulationResult
  trajectories::Dict{ASCIIString,PopulationTrajectory}
  metadata::Dict{ASCIIString,Any}
end

function SimulationResult(x::Vector{Species})
  d = Dict{ASCIIString,PopulationTrajectory}()
  for i in eachindex(x)
    d[x.id] = PopulationTrajectory(x)
  end
  return SimulationResult(d, Dict{ASCIIString,Any}())
end

function SimulationResult(x::Vector{Species}, n::Int)
  d = Dict{ASCIIString,PopulationTrajectory}()
  for i in eachindex(x)
    d[x.id] = PoulationTrajectory(x.id, n)
  end
  return SimulationResult(d, Dict{ASCIIString,Any}())
end

function getindex(sr::SimulationResult, key...)
  getindex(srt.trajectories, value, key)
end

function setindex!(sr::SimulationResult, value, key...)
  setindex!(sr.trajectories, value, key)
end
