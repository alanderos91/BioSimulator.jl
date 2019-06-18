struct SamplePath{T1,T2}
  t::T1
  x::T2
end

Ensemble{T1,T2} = Vector{SamplePath{T1,T2}}

SamplePath{T1,T2}() where {T1, T2}  = SamplePath(T1(), T2())

# SamplePath() = SamplePath{Float64,Int}()

# Ensemble(ntrials) = [SamplePath() for i in 1:ntrials]

function Base.push!(xw::SamplePath, t, x)
  push!(xw.t, t)

  if eltype(xw.x) <: AbstractVector
    push!(xw.x, copy(x))
  elseif eltype(xw.x) <: Tuple
    push!(xw.x, tuple(x...))
  end

  return xw
end

function Base.sizehint!(xw::SamplePath, n)
  sizehint!(xw.t, n)
  sizehint!(xw.x, n)

  return xw
end

function update!(xw::SamplePath, t, x)
  @inbounds push!(xw, t, x)
end
