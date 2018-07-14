struct AveragePath{T}
  tdata  :: Vector{T}
  xmean  :: Matrix{T}
  xstd   :: Matrix{T}
end

function AveragePath(xw :: RegularEnsemble{T1,T2}) where {T1,T2}
  tdata = xw[1].tdata

  d = size(xw[1].xdata, 1)
  n = length(tdata)
  trials = length(xw)

  xmean = mean(xw[k].xdata for k in 1:trials)
  temp = std([xw[k].xdata[i, j] for i in 1:d, j in 1:n, k in 1:trials], 3)
  xstd = reshape(temp, n, d)

  return AveragePath{T1}(tdata, xmean, xstd)
end

@recipe function f(xw :: AveragePath)
  tdata  = xw.tdata
  xmean  = xw.xmean'
  xstd   = xw.xstd

  legend     --> true
  grid       --> false
  xlims      --> (tdata[1], tdata[end])
  ylims      --> (0.0, Inf)
  seriestype --> :scatter
  yerrorbar  --> xstd

  tdata, xmean
end
