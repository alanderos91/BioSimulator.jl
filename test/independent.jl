T = 10.
x0 = 1_000
model_size = [10, 100, 500]

ssa = SSA(end_time=T)
frm = FRM(end_time=T)
nrm = NRM(end_time=T)
odm = ODM(end_time=T)
sal = SAL(end_time=T)

algorithms = [ssa, frm, nrm, odm, sal]

for algorithm in algorithms
    @printf "%+6s\n" uppercase(string(typeof(algorithm)))
    for n in model_size
        model = independent(n, x0)

        @printf "%+6s: %3d" "n" n
        @time simulate(model, algorithm)
    end
end
