T = 10.
x0 = 1_000
model_size = [10, 100, 500]

ssa = SSA(T)
frm = FRM(T)
nrm = NRM(T)
odm = ODM(T, 1000)
sal = SAL(T, 0.125, 100.0, 0.75)
algorithms = [ssa, frm, nrm, odm, sal]

for algorithm in algorithms
    @printf "%+6s\n" uppercase(string(typeof(algorithm)))
    for n in model_size
        model = independent(n, x0)

        @printf "%+6s: %3d" "n" n
        @time simulate(model, algorithm)
    end
end
