using BioSimulator
include("test_models.jl")

function test()
    for alg in BioSimulator.ALGORITHMS
        @time simulate(kendall(), with=alg)
    end
    Profile.clear_malloc_data()
    for alg in BioSimulator.ALGORITHMS
        @time simulate(kendall(), T=4.0, dt=0.5, with=alg, itr=100_000)
    end
end

test()
