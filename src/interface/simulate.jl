function simulate(algorithm, X0, r; time::Float64=1.0, sampling_interval::Float64=1.0, nrlz::Integer=1)

  Xt = similar(X0)

  npts = round(Int, time / sampling_interval + 1)

  output = PartialHistory(length(Xt), npts, nrlz, 0.0, time)

  simulate!(output, Xt, algorithm, X0, r, nrlz)
end

function simulate!(output, Xt, algorithm, X0, r, nrlz) # nrlz is encoded in PartialHistory; refactor

  for i in 1:nrlz
    init!(algorithm)
    copy!(Xt, X0)
    compute_propensities!(r, Xt)
    interval = 1

    while !done(algorithm)
      output, interval = update!(output, Xt, get_time(algorithm), interval, i)
      step!(algorithm, Xt, r)
    end

    output, interval = update!(output, Xt, get_time(algorithm), interval, i)
  end

  return output
end
