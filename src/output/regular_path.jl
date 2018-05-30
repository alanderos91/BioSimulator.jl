struct RegularPath{T1,T2}
  tdata :: Vector{T1}
  xdata :: Vector{Vector{T2}}
end

RegularEnsemble{T1,T2} = Vector{RegularPath{T1,T2}}

# constructors

@inline function RegularPath{T1,T2}(nspecies::Int, epochs::Int) where {T1,T2}
  n = epochs + 1
  tdata = zeros(T1, n)
  xdata = [zeros(T2, nspecies) for i in 1:n]

  return RegularPath(tdata, xdata)
end

RegularPath(nspecies, epochs) = RegularPath{Float64,Int}(nspecies, epochs)

RegularEnsemble(ntrials, nspecies, epochs) = [RegularPath(nspecies, epochs) for i in 1:ntrials]

function update!(xw :: RegularPath, t, x, epoch)
  tdata = xw.tdata
  xdata = xw.xdata

  max_epoch = length(tdata)

  while (epoch <= max_epoch) && (t >= tdata[epoch])
    @inbounds copy!(xdata[epoch], x)
    epoch = epoch + 1
  end

  return epoch
end
