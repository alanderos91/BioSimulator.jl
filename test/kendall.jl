using BioSimulator
using Base.Test

function kendall_mean(i,t,alpha,mu,nu)
   x = exp((alpha-mu)*t)
  return i * x + nu/(alpha-mu)*(x-1)
end

function computed_mean(sj)
  mval = 0.0
  for r in sj.results
    mval = mval + r["Particle"].states[end].pop
  end
  return mval / length(sj.results)
end

x = [Species("Particle", 5, true)]
r = [Reaction("Birth", "alpha", [1], [2]),
     Reaction("Death", "mu", [1], [0]),
     Reaction("Immigration", "nu", [0], [1])]
p = Dict{ASCIIString, Float64}("alpha" => 2.0,
                               "mu" => 1.0,
                               "nu" => 0.5)

kendall = Network("Kendall's Process", x, r, p);

t_final = 4.0
Δt = 0.1
seed = 5357

# Compute mean for comparison
points = round(Int, t_final/Δt) + 1
t = linspace(0.0,t_final, points)
theoretical = kendall_mean(x[1].pop,t,p["alpha"],p["mu"],p["nu"])

algorithms = [:ssa, :odm, :nrm, :sal, :frm]

# Run SSA and SAL once to compile
print("Precompiling..."); @time begin
  for a in algorithms
    simulate(Simulation(kendall), t_final, a, o=Explicit(), itr=1)
    simulate(Simulation(kendall), t_final, a, o=Uniform(), dt=0.1, itr=1)
  end
end

println("Running tests...")
for a in algorithms
  print(" * Explicit ", uppercase(string(a)))
  srand(seed); @time result = simulate(Simulation(kendall), t_final, a, o=Explicit(), itr=10^5)
  computed = computed_mean(result)
  @test_approx_eq_eps computed theoretical[end] 1e0

  print(" *  Uniform ", uppercase(string(a)))
  srand(seed); @time result = simulate(Simulation(kendall), t_final, a, o=Uniform(), dt=0.1, itr=10^5)
  computed = computed_mean(result)
  @test_approx_eq_eps computed theoretical[end] 1e0
end
