import NewBioSimulator: parse_model

@testset "birth-death-immigration" begin
  ALGORITHMS = [
    (Direct(), HasRates),
    (Direct(), HasSums),
    (EnhancedDirect(), HasRates),
    (SortingDirect(), HasRates),
    (FirstReaction(), HasRates),
    (NextReaction(), HasRates),
    (RejectionSSA(), HasRates)
  ]

  N = 10_000

  function kendall_mean(i, t, α, μ, ν)
    x = exp((α - μ) * t)
    return i * x + ν / (α - μ) * (x - 1)
  end

  state, model = parse_model(kendall())

  expected = kendall_mean(5, 4.0, 2.0, 1.0, 0.5)

  @testset "$(alg), $(rates_cache)" for (alg, rates_cache) in ALGORITHMS
    msg = rates_cache == HasRates ? "linear search" : "binary search"

    @info "Precompiling $(alg) using $(msg)...\n"
    @time simulate(state, model, alg, 4.0, rates_cache)

    @info "Running $(alg) using $(msg)...\n"
    @time result = [simulate(state, model, alg, 4.0, rates_cache)[end][1] for i in 1:N]

    println("  absolute error = $(abs(mean(result) - expected))\n")
  end
end
