immutable SimulationResult
  trajectories::Dict{ASCIIString,PopulationTrajectory}
  metadata::Dict{ASCIIString,Any}
end

function SimulationResult(x::Vector{Species})
  d = Dict{ASCIIString,PopulationTrajectory}()
  for i in eachindex(x)
    d[x[i].id] = PopulationTrajectory(x[i])
  end
  return SimulationResult(d, Dict{ASCIIString,Any}())
end

function SimulationResult(x::Vector{Species}, n::Int)
  d = Dict{ASCIIString,PopulationTrajectory}()
  for i in eachindex(x)
    d[x[i].id] = PopulationTrajectory(x[i].id, n)
  end
  return SimulationResult(d, Dict{ASCIIString,Any}())
end

function getindex(sr::SimulationResult, key)
  sr.trajectories[key]
end

function setindex!(sr::SimulationResult, value, key)
  sr.trajectories[key] = value
end

function update!(sr::SimulationResult, t::Float64, spcs::Vector{Species})
  for s in spcs
    update!(sr[s.id], t, s)
  end
end

function update!(sr::SimulationResult, t::Float64, spcs::Vector{Species}, j::Int)
  for s in spcs
    update!(sr[s.id], t, s, j)
  end
end
