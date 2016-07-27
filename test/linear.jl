T = 10.
x0 = 1_000
model_size = [10, 100, 500]

for a in BioSimulator.ALGORITHMS
    @printf "%+6s\n" uppercase(string(a))
    for M in model_size
        model = linear(M, x0)
        @printf "%+6s: %3d" "M" M
        @time simulate(model, time=T, method=a, track=[:S1])
    end
end
