type PartialHistory
  t     :: LinSpace{Float64}
  data  :: Array{Int, 3}
end

function PartialHistory(
    d1      :: Integer, # number of species
    d2      :: Integer, # number of intervals
    d3      :: Integer, # number of realizations
    start   :: Float64, # start time
    stop    :: Float64 # end time
)

  return PartialHistory(linspace(start, stop, d2), zeros(Int, d1, d2, d3))
end

function update!(
    x           :: PartialHistory,
    Xt          :: Vector,
    t           :: AbstractFloat,
    interval    :: Integer,
    realization :: Integer
)

  ttt      = x.t
  data     = x.data
  maxindex = length(ttt)

  while (interval <= maxindex) && (t >= ttt[interval])
    for k in eachindex(Xt)
      data[k, interval, realization] = Xt[k]
    end
    interval = interval + 1
  end

  return interval
end
