struct SamplePath{T1,T2}
  tdata :: Vector{T1}
  xdata :: Vector{Vector{T2}}
end

Ensemble{T1,T2} = Vector{SamplePath{T1,T2}}

# constructors

SamplePath{T1,T2}() where {T1, T2} = SamplePath(T1[], Vector{Vector{T2}}())
SamplePath() = SamplePath{Float64,Int}()

Ensemble(ntrials) = [SamplePath() for i in 1:ntrials]

function Base.push!(xw :: SamplePath, t, x)
  push!(xw.tdata, t)
  push!(xw.xdata, copy(x))

  return xw
end

function Base.sizehint!(xw :: SamplePath, n)
  sizehint!(xw.tdata, n)
  sizehint!(xw.xdata, n)

  return xw
end

function update!(xw :: SamplePath, t, x)
  @inbounds push!(xw, t, x)
end
