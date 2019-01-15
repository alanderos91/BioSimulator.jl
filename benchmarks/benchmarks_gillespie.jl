using BenchmarkTools, Gillespie, StatsBase

function ssa_wrapper(x, x0, F, nu, parms, t_final, n_trial)
  for i in 1:n_trial
    copy!(x, x0)
    ssa(x, F, nu, parms, t_final)
  end

  return nothing
end

function summarize(results, benchmark_set)

  for m in benchmark_set
    file = open("./gillespie/$(m).txt", "a")

    result = results[m]

    println(file, "model: $(m)")
    println(file)
    println(file, "n = $(length(result.times))")
    println(file, "mean:               $(mean(result.times) * 1e-6) ms")
    println(file, "std:                $(std(result.times) * 1e-6) ms")
    println(file, "median:             $(median(result.times) * 1e-6) ms")
    println(file, "mad (normal=true):  $(mad(result.times, normalize=true) * 1e-6) ms")
    println(file, "mad (normal=false): $(mad(result.times, normalize=false) * 1e-6) ms")
    println(file)

    close(file)
  end

  return nothing
end

# initialize test suite
suite = BenchmarkGroup()

# load the following variables:
# - parameters
# - MODELS
# - SEED
include("benchmark_parameters.jl")

# build handy collections
const benchmark_set = intersect(ARGS, MODELS)

# build up the benchmark suite
for m in benchmark_set
  # unpack parameters
  t_final, n_saves, n_trial, n_sample, t_limit = parameters[m]
      
  # load model
  include(joinpath(m, "gillespie.jl"))

  x = similar(x0)

  suite[m] = @benchmarkable ssa_wrapper($x, $x0, $F, $nu, $parms, $t_final, $n_trial) samples=n_sample seconds=t_limit setup=(srand($SEED))
end

results = run(suite, verbose = true)

# summarize(results, benchmark_set)
BenchmarkTools.save("gillespie.json", results)
