using BioSimulator

include("test_models.jl")

function test(algorithms)

    for alg in algorithms
        @time simulate(kendall(), with=alg)
    end

    Profile.clear_malloc_data()

    for alg in algorithms
        @time simulate(kendall(), time=4.0, sampling_interval=0.5, method=alg, realizations=100_000)
    end

    nothing
end

cmdargs = map(symbol, ARGS)

test(intersect(BioSimulator.ALGORITHMS, cmdargs))
