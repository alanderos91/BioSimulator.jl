struct AveragePath{T}
  tdata  :: Vector{T}
  xmean  :: Matrix{T}
  xstd   :: Matrix{T}
end

function AveragePath(xw :: RegularEnsemble{T1,T2}) where {T1,T2}
  tdata = xw[1].tdata
  trials = length(xw)

  # temporary fix
  x = [xw[k].xdata for k in 1:trials]

  return AveragePath{T1}(tdata, mean(x), std(x))
end

@recipe function f(xw :: AveragePath)
  tdata  = xw.tdata
  xmean  = transpose(xw.xmean)
  xstd   = transpose(xw.xstd)

  legend     --> true
  grid       --> false
  xlims      --> (tdata[1], tdata[end])
  ylims      --> (0.0, Inf)
  seriestype --> :scatter
  yerrorbar  --> xstd

  tdata, xmean
end
