using BenchmarkTools, BioSimulator, StatsBase

function summarize(results, benchmark_set)

  for m in benchmark_set
    fname = Threads.nthreads() > 1 ? "biosimulator-parallel/$(m).txt" : "biosimulator-serial/$(m).txt"
    
    file = open(fname, "a")

    for a in ALGOS
      for o in OPTIONS
        result = results[m][a][o]

        println(file, "model: $(m) / algorithm: $(a) / output: $(o)")
        println(file)
        println(file, "n = $(length(result.times))")
        println(file, "mean:               $(mean(result.times) * 1e-6) ms")
        println(file, "std:                $(std(result.times) * 1e-6) ms")
        println(file, "median:             $(median(result.times) * 1e-6) ms")
        println(file, "mad (normal=true):  $(mad(result.times, normalize=true) * 1e-6) ms")
        println(file, "mad (normal=false): $(mad(result.times, normalize=false) * 1e-6) ms")
        println(file)
      end
    end
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
const algtypemap = Dict(
                  "dm"  => Direct(),
                  "frm" => FirstReaction(),
                  "nrm" => NextReaction(),
                  "odm" => OptimizedDirect(),
                  "otl" => TauLeaping(),
                  "sal" => StepAnticipation()
                )

const outtypemap = Dict("fixed" => Val(:fixed), "full" => Val(:full))

const ALGOS   = collect(keys(algtypemap))
const OPTIONS = collect(keys(outtypemap))
const benchmark_set = intersect(ARGS, MODELS)

# build up the benchmark suite
for m in benchmark_set
  # build the set for the model m
  suite[m] = BenchmarkGroup()

  for a in ALGOS
    # build the set for the algorithm a
    suite[m][a] = BenchmarkGroup()

    # select algorithm
    algtype = algtypemap[a]

    for o in OPTIONS
      # select output option
      outtype = outtypemap[o]

      # unpack parameters
      t_final, n_saves, n_trial, n_sample, t_limit = parameters[m]
      
      # load model
      include(joinpath(m, "biosimulator.jl"))

      suite[m][a][o] =
        @benchmarkable simulate($model, $algtype, $outtype,
          time   = $t_final,
          epochs = $n_saves,
          trials = $n_trial) samples=n_sample seconds=t_limit setup=(srand($SEED))
    end
  end
end

results = run(suite, verbose = true)

# summarize(results, benchmark_set)

fname = Threads.nthreads() > 1 ? "biosimulator-parallel.json" : "biosimulator-serial.json"

BenchmarkTools.save(fname, results)
