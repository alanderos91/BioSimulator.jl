using BioSimulator

include("test_models.jl")

function test(algorithms)

    for alg in algorithms
        @time simulate(autoreg(), alg, time=1.0, epochs=1, trials=1)
    end

    Profile.clear_malloc_data()

    for alg in algorithms
        @time simulate(autoreg(), alg, time=1000.0, epochs=1000, trials=10)
    end

    nothing
end

cmdargs = map(Symbol, ARGS)
algorithms = []

if :SSA ∈ cmdargs push!(algorithms, SSA) end
if :FRM ∈ cmdargs push!(algorithms, FRM) end
if :ODM ∈ cmdargs push!(algorithms, ODM) end
if :NRM ∈ cmdargs push!(algorithms, NRM) end
if :SAL ∈ cmdargs push!(algorithms, SAL) end

test(algorithms)
