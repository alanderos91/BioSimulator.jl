function kendall_mean(i, t, α, μ, ν)
  x = @. exp((α - μ) * t)
  @. x = i * x + ν / (α - μ) * (x - 1)
  return x
end

i = 5
α = 2.0
μ = 1.0
ν = 0.5

t = 4.0
n = 100
u = 5357
m = 1_000 # for a real test, use 100_000

model = kendall(i, α, μ, ν)

theoretical = kendall_mean(i, range(0.0, stop = t, length = n + 1), α, μ, ν)

algorithms = [
  Direct(),
  FirstReaction(),
  NextReaction(),
  OptimizedDirect(),
  TauLeaping(),
  StepAnticipation()
]

# Run SSA and SAL once to compile
print("    Precompiling..."); @time begin
  for algorithm in algorithms
    run_test(model, algorithm, t, n, 1)
  end
end

print("    Running tests...\n\n")
for algorithm in algorithms
  print("   - $(algorithm): ")

  seed!(u)

  @time result = run_test(model, algorithm, t, n, m)

  # count the number of relative errors that lie outside the interval [0.98, 1.02]
  avg = BioSimulator.AveragePath(result.simulation_data)
  
  observed = reshape(avg.xmean, n+1, 1)
  relative = observed ./ theoretical
  badness = count(x -> !isapprox(x, 1.0, rtol=0.4), relative)

  @test badness / n ≤ 0.05

  print("     % bad estimates = ", badness / m, " out of $(m)\n")
  println()
end
