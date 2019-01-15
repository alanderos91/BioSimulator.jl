using BenchmarkTools, StatsBase

##### function definitions #####

# for BenchmarkTools.jl results
function summarize(::Val{true}, fname)
  # retrieve the Trial objects for the test suite
  results = BenchmarkTools.load(fname)[1]

  for (b_info, trial) in leaves(results)
    # convert results from ns to ms
    times = [ t * 1e-6 for t in trial.times ]

    b_mean   = mean(times)
    b_std    = std(times)
    b_median = median(times)
    b_q1     = quantile(times, 0.25)
    b_q3     = quantile(times, 0.75)

    # benchmark info string
    str = [ "$(s) " for s in b_info ]

    println("----- $(str...)-----")
    println("trials: $(length(times))")
    println("  mean: $(b_mean) ms")
    println("   std: $(b_std) ms")
    println("median: $(b_median) ms")
    println("    Q1: $(b_q1) ms")
    println("    Q3: $(b_q3) ms")
    println()
  end
end

function summarize(::Val{false}, fname)
  # retrieve data from Python benchmarks
  results = map(x -> parse(Float64, x), readlines(fname))

  # convert results from s to ms
  times = [ t * 1e3 for t in results ]

  b_mean   = mean(times)
  b_std    = std(times)
  b_median = median(times)
  b_q1     = quantile(times, 0.25)
  b_q3     = quantile(times, 0.75)

  println("----- $(fname) ------")
  println("trials: $(length(results))")
  println("  mean: $(b_mean) ms")
  println("   std: $(b_std) ms")
  println("median: $(b_median) ms")
  println("    Q1: $(b_q1) ms")
  println("    Q3: $(b_q3) ms")
  println()
end

##### command-line tool #####

fname = ARGS[1]
ftype = split(fname, ".")[end]

is_json = ftype == "json"

summarize(Val(is_json), fname)
