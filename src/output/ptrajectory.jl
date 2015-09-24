import Base: length, show, getindex, push!, setindex!, isempty

immutable PopulationTrajectory
  id::ASCIIString
  states::Vector{PopulationState}
end

function PopulationTrajectory(id::ASCIIString)
  return PopulationTrajectory(id, PopulationState[])
end

function PopulationTrajectory(id::ASCIIString, n::Int)
  n < 0 && error("invalid PopulationTrajectory length")
  return PopulationTrajectory(id, Array{PopulationState}(n))
end

function PopulationTrajectory(x::Species; t_start::Float64=0.0)
  state = PopulationState(t_start, x.pop)
  return PopulationTrajectory(x.id, [state])
end

push!(tr::PopulationTrajectory, s::PopulationState) = push!(tr.states, s)
getindex(tr::PopulationTrajectory, i::Integer) = getindex(tr.states, i)

function setindex!(tr::PopulationTrajectory, s::PopulationState, i::Integer)
  setindex!(tr.states, s, i)
end

function isempty(tr::PopulationTrajectory)
  return isempty(tr.states)
end

function length(tr::PopulationTrajectory)
  return length(tr.states)
end

function get_id(tr::PopulationTrajectory)
  return tr.id
end

function show(io::IO, tr::PopulationTrajectory)
  @printf io "   Name        Time       Value  \n"
  @printf io "----------  ----------  ---------\n"
  for state in t.states
    @printf " %8s   " tr.id
    show(io, state)
  end
end

function update!(tr::PopulationTrajectory, t::Float64, s::Species)
  push!(tr, PopulationState(t, s.pop))
end

function discretize(tr::PopulationTrajectory; dt::Float64=1.0, t_start::Float64=0.0)
  t_final = tr[end].tval
  n = round(Int, (t_final - t_start) / dt) + 1
  ntr = PopulationTrajectory(tr.id, n)
  t_next = t_start
  j = 1

  for state in tr.states
    while state.tval >= t_next
      ntr[j] = PopulationState(tr.tval, tr.pop)
      j = j + 1
      t_next = t_next + dt
    end
  end
  return ntr
end
