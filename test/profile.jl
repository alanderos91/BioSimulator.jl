using BioSimulator

include("test_models.jl")

function test(algorithms)

    for alg in algorithms
        @time simulate(autoreg(), method=alg)
    end

    Profile.clear_malloc_data()

    for alg in algorithms
        @time simulate(autoreg(), time=1000.0, sampling_interval=10.0, method=alg, realizations=100)
    end

    nothing
end

cmdargs = map(symbol, ARGS)

test(intersect(BioSimulator.ALGORITHMS, cmdargs))
