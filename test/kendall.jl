using BioSimulator
using Base.Test
using DataFrames

function kendall_mean(i,t,α,μ,ν)
   x = exp((α-μ)*t)
  return i * x + ν/(α-μ)*(x-1)
end

m = kendall()
x = m.species[:X].population
α = p[:α].value
μ = p[:μ].value
ν = p[:ν].value
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
    simulate(m, T=T, with=a, output=Explicit(), itr=1)
    simulate(m, T=T, with=a, output=Uniform(), dt=dt, itr=1)
  end
end

print("    Running tests...\n\n")
for a in BioSimulator.ALGORITHMS
  # print(" * Explicit ", uppercase(string(a)))
  # srand(seed); @time result = simulate(Simulation(m), T, a, o=Explicit(), itr=100_000)
  # computed = computed_mean(result)
  # @test_approx_eq_eps computed theoretical[end] 1e0

  print("   - Uniform ", uppercase(string(a)))
  srand(seed); @time result = simulate(m, T=T, with=a, output=Uniform(), dt=dt, itr=itr)
  computed = aggregate(species(result), :time, mean)
  print("     |observed - theoretical| = ", abs(computed[:X_mean][end] - theoretical[end]), "\n")
  println()
end
