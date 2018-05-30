x0 = 1_000
model_size = [10, 100, 500]

t = 10.0
n = 1
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
  @printf "%+6s\n" algorithm
  for M in model_size
      model = independent(M, x0)
      srand(u)
      @printf "%+6s: %3d" "M" M
      @time run_test(model, algorithm, t, n, m)
  end
end
