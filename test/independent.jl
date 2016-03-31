T = 10
x0 = 1_000
model_size = [10, 100, 500]

for a in BioSimulator.ALGORITHMS
    @printf "%+6s\n" uppercase(string(a))
    for n in model_size
        model = independent(n, x0)
        @printf "%+6s: %d" "n" n
        @time simulate(model, T=T, with=a, track=[:S1])
    end
end
