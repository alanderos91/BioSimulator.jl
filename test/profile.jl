using BioSimulator

include("test_models.jl")

function test(algorithms)

    for alg in algorithms
        @time simulate(autoreg(), alg)
    end

    Profile.clear_malloc_data()

    for alg in algorithms
        @time simulate(autoreg(), alg, sampling_interval=10.0, nrlz=50)
    end

    nothing
end

cmdargs = map(symbol, ARGS)
algorithms = []

if :SSA ∈ cmdargs push!(algorithms, SSA(1000.0)) end
if :FRM ∈ cmdargs push!(algorithms, FRM(1000.0)) end
if :ODM ∈ cmdargs push!(algorithms, ODM(1000.0, 1000)) end
if :NRM ∈ cmdargs push!(algorithms, NRM(1000.0)) end
if :SAL ∈ cmdargs push!(algorithms, SAL(1000.0, 0.125, 100.0, 0.75)) end

test(algorithms)
