
# built from RecursiveArrayTools.jl
# based on DiffEqArray struct
struct SamplePath{T,N,A,B} <: AbstractVectorOfArray{T, N}
  u::A # A <: AbstractVector{<: AbstractArray{T, N - 1}}
  t::B
end

Ensemble{T,N,A,B} = Vector{SamplePath{T,N,A,B}}

SamplePath(xs::AbstractVector{T}, ts, dims::NTuple{N}) where {T, N} = SamplePath{eltype(T), N, typeof(xs), typeof(ts)}(xs, ts)
# Assume that the first element is representative all all other elements
SamplePath(xs::AbstractVector, ts::AbstractVector) = SamplePath(xs, ts, (size(xs[1])..., length(xs)))

@inline Base.getindex(xw::SamplePath{T, N}, I::AbstractArray{Int}) where {T, N} = SamplePath(xw.u[I], xw.t[I])

Base.show(io::IO, xw::SamplePath) = (print(io,"t: ");show(io, xw.t);println(io);print(io,"x: ");show(io, xw.u))
Base.show(io::IO, m::MIME"text/plain", xw::SamplePath) = (print(io,"t: ");show(io,m,xw.t);println(io);print(io,"x: ");show(io,m,xw.u))

# SamplePath() = SamplePath{Float64,Int}()

# Ensemble(ntrials) = [SamplePath() for i in 1:ntrials]

function Base.sizehint!(xw::SamplePath, n)
  sizehint!(xw.t, n)
  sizehint!(xw.u, n)

  return xw
end

function update!(xw::SamplePath, t, x, save_points)
  i = searchsortedlast(save_points, t)

  if !(save_points[i] in xw.t)
    push!(xw.u, copy(x))        # update the data
    push!(xw.t, save_points[i]) # update the time series
  end

  return xw
end

function update!(xw::SamplePath, t, x, save_points::Nothing)
  push!(xw.u, copy(x)) # update the data
  push!(xw.t, t)       # update the time series

  return xw
end

##### interpolation
function get_regular_path(xw::SamplePath, tfinal, epochs) where {T1,T2}
  max_epoch = epochs + 1

  if xw.u[1] isa Number
    x0 = [xw.u[1]]
  else
    x0 = xw.u[1]
  end

  ts = collect(range(0.0, stop = tfinal, length = max_epoch))
  new_xw = SamplePath([x0], ts)
  epoch = 2

  # copy data from the sample path
  for j in 2:length(xw)
    @inbounds t = xw.t[j]
    @inbounds x = xw.u[j]

    while (epoch <= max_epoch) && (t >= ts[epoch])
      push!(new_xw.u, x)
      epoch = epoch + 1
    end
  end

  # fill in the regular path
  while epoch <= max_epoch
    x = xw.u[end]
    push!(new_xw.u, x)
    epoch = epoch + 1
  end

  return new_xw
end

function get_regular_ensemble(x::Ensemble, tfinal, epochs)
  return [get_regular_path(x[i], tfinal, epochs) for i in eachindex(x)]
end

##### plotting recipes #####

@recipe function f(xw::SamplePath)
  seriestype --> :steppre

  xw.t, xw'
end

@recipe function f(xw::SamplePath{T,1}) where {T}
  seriestype --> :steppre

  xw.t, xw.u
end

@recipe function f(ens::Ensemble, epochs = 100)
  # regularize the sample paths
  tfinal = ens[1].t[end]

  reg = get_regular_ensemble(ens, tfinal, epochs)

  # extract the series data
  ts = reg[1].t
  xs = convert(Array, mean(reg)')
  bars = convert(Array, std(reg)')

  # make the default a scatter plot
  # and add bars to represent standard deviation
  seriestype --> :scatter
  yerrorbar  --> bars

  ts, xs
end
