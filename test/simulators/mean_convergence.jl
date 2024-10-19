import BioSimulator: parse_model

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

  x0 = 5
  birth_rate = 2.0
  death_rate = 1.0
  immgr_rate = 0.5
  tfinal = 4.0
  delta_t = tfinal

  state, model = parse_model(kendall(
    x0 = x0,
    birth_rate = birth_rate,
    death_rate = death_rate,
    immgr_rate = immgr_rate)
  )

  expected = kendall_mean(x0, tfinal, birth_rate, death_rate, immgr_rate)
  save_points = range(0, tfinal, step = delta_t)
  @testset "$(alg), $(rates_cache)" for (alg, rates_cache) in ALGORITHMS
    msg = (rates_cache == HasRates) ? "linear search" : "binary search"

    @info "Precompiling $(alg) using $(msg)...\n"
    @time simulate(state, model, alg, tfinal = tfinal, save_points = save_points, ntrials = 1, rates_cache = rates_cache)

    @info "Running $(alg) using $(msg)...\n"
    @time result = simulate(state, model, alg, tfinal = tfinal, save_points = save_points, ntrials = N, rates_cache = rates_cache)

    println("  absolute error = $(abs(mean(result).u[end][1] - expected))\n")
  end

  TAULEAPING = [
    TauLeapingDG2001(), TauLeapingDGLP2003(),
    StepAnticipation(), HybridSAL()
  ]

  @testset "$(alg)" for alg in TAULEAPING
    @info "Precompiling $(alg)...\n"
    @time simulate(state, model, alg, tfinal = tfinal, save_points = save_points, ntrials = 1)
    @info "Running $(alg)...\n"
    @time result = simulate(state, model, alg, tfinal = tfinal, save_points = save_points, ntrials = N)

    println("  absolute error = $(abs(mean(result).u[end][1] - expected))\n")
  end
end
