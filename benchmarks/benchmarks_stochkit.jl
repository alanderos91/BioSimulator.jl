using BenchmarkTools, StatsBase

function summarize(results, benchmark_set)

  for m in benchmark_set
    file = open("./stochkit/$(m).txt", "a")

    result = results[m]
    adjustment = adjust[m]

    for i in eachindex(result.times)
      # convert to ms and correct for overhead
      result.times[i] = result.times[i] * 1e-6 - adjustment
    end

    println(file, "model: $(m)")
    println(file)
    println(file, "n = $(length(result.times))")
    println(file, "mean:               $(mean(result.times)) ms")
    println(file, "std:                $(std(result.times)) ms")
    println(file, "median:             $(median(result.times)) ms")
    println(file, "mad (normal=true):  $(mad(result.times, normalize=true)) ms")
    println(file, "mad (normal=false): $(mad(result.times, normalize=false)) ms")
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
const adjust = Dict(
  "kendall"     => 80.0,
  "mmek"        => 100.0,
  "autoreg"     => 115.0,
  "dimer-decay" => 100.0,
  "yeast"       => 130.0
)

# build up the benchmark suite
for m in benchmark_set
  # unpack parameters
  t_final, n_saves, n_trial, n_sample, t_limit = parameters[m]
  
  # build filename
  fname = "./$(m)/stochkit2.xml"

  # build command
  cmd = `ssa -m $(fname) -t $(t_final) -i $(n_saves) -r $(n_trial) -f --seed $(SEED) '&>' /dev/null`

  suite[m] = @benchmarkable run($cmd) samples=n_sample seconds=t_limit
end

results = run(suite, verbose = true)

# perform adjustments
for b in leaves(results)
  b_info, trial = b

  m = b_info[1]

  adjustment = adjust[m]

  for i in eachindex(trial.times)
    trial.times[i] -= adjustment
  end
end

# summarize(results, benchmark_set)
BenchmarkTools.save("stochkit.json", results)
