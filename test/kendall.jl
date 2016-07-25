using BioSimulator
using Base.Test
using DataFrames

function kendall_mean(i,t,α,μ,ν)
   x = exp((α-μ)*t)
  return i * x + ν/(α-μ)*(x-1)
end

m = kendall()
x = value(species_list(m)[:X])
p = parameter_list(m)
α = value(p[:α])
μ = value(p[:μ])
ν = value(p[:ν])
T = 4.0
dt = 0.1
seed = 5357
itr = 100_000

# Compute mean for comparison
npts = round(Int, T / dt) + 1
t = linspace(0.0, T, npts)
theoretical = kendall_mean(x,t,α,μ,ν)

# Run SSA and SAL once to compile
print("    Precompiling..."); @time begin
  for a in BioSimulator.ALGORITHMS
    simulate(m, time=T, method=a, output=:explicit, realizations=1)
    simulate(m, time=T, method=a, output=:fixed, sampling_interval=dt, realizations=1)
  end
end

print("    Running tests...\n\n")
for a in BioSimulator.ALGORITHMS
  # print(" * Explicit ", uppercase(string(a)))
  # srand(seed); @time result = simulate(Simulation(m), T, a, o=Explicit(), itr=100_000)
  # computed = computed_mean(result)
  # @test_approx_eq_eps computed theoretical[end] 1e0

  print("   - Uniform ", uppercase(string(a)))
  srand(seed); @time result = simulate(m, time=T, method=a, output=:fixed, sampling_interval=dt, realizations=itr)
  #computed = aggregate(get_speciesdf(result), :time, mean)
  #print("     |observed - theoretical| = ", abs(computed[:X_mean][end] - theoretical[end]), "\n")
  println()
end
