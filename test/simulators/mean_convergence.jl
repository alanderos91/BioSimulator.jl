import NewBioSimulator: parse_model

@testset "birth-death-immigration" begin
  N = 100_000

  function kendall_mean(i, t, α, μ, ν)
    x = exp((α - μ) * t)
    return i * x + ν / (α - μ) * (x - 1)
  end

  state, model = parse_model(kendall())

  expected = kendall_mean(5, 4.0, 2.0, 1.0, 0.5)
  result   = zeros(Int, N)

  @testset "$(alg), $(rates_cache)" for (alg, rates_cache) in [(Direct(), HasRates), (Direct(), HasSums), (FirstReaction(), HasRates)]
    msg = rates_cache == HasRates ? "linear search" : "binary search"

    @info "Precompiling $(alg) using $(msg)...\n"
    @time simulate(state, model, alg, 4.0, rates_cache)

    @info "Running $(alg) using $(msg)...\n"
    @time for i in 1:N
      result[i] = simulate(state, model, alg, 4.0, rates_cache)[1]
    end
    # @test mean(result) ≈ expected
    println("  absolute error = $(abs(mean(result) - expected))")
  end
end