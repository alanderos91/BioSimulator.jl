struct BioSimHistogram{T}
  xdata :: Matrix{T}
end

function BioSimHistogram(xw :: RegularEnsemble{T1,T2}) where {T1,T2}
  xdata = transpose(hcat( (xw[i].xdata[:, end] for i in 1:length(xw))...))

  return BioSimHistogram{T2}(xdata)
end

function BioSimHistogram(xw :: Ensemble{T1,T2}) where {T1,T2}
  xdata = transpose(hcat( (xw[i].xdata[end] for i in 1:length(xw))... ))

  return BioSimHistogram{T2}(xdata)
end

@recipe function f(xw :: BioSimHistogram)
  xdata = xw.xdata

  legend     --> true
  seriestype --> :histogram

  xdata
end
