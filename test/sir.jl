#=
Updating propensities when S population goes to 0 may lead to negative propensities. This happened when we failed to properly update the propensities using Kahan summation. This test only checks that we do not throw an error, but it should always run *except* if the propensities were not updated correctly.
=#
S = 9999
I = 1
R = 0

model = sir(S,I,R)

t = 1000.0
n = 1000
u = 5357
m = 1

algorithms = [
  Direct(),
  FirstReaction(),
  NextReaction(),
  OptimizedDirect(),
  TauLeaping(),
  StepAnticipation()
]

for algorithm in algorithms
  seed!(u)
    print("   - $(algorithm): ")
  @time run_test(model, algorithm, t, n, m)
end
