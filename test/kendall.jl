using BioSimulator
using Base.Test
using DataFrames

function kendall_mean(i,t,alpha,mu,nu)
   x = exp((alpha-mu)*t)
  return i * x + nu/(alpha-mu)*(x-1)
end

m = kendall()
x = m.species[:X].population
t_final = 4.0
Δt      = 0.1
seed    = 53571

# Compute mean for comparison
points = round(Int, t_final/Δt) + 1
t      = linspace(0.0,t_final, points)
theoretical = kendall_mean(x,t,α,μ,ν)

algorithms = [:ssa, :odm, :nrm, :sal, :frm]

# Run SSA and SAL once to compile
print("    Precompiling..."); @time begin
  for a in algorithms
    simulate(m, tf=t_final, with=a, output=Explicit(), itr=1)
    simulate(m, tf=t_final, with=a, output=Uniform(), dt=0.1, itr=1)
  end
end

print("    Running tests...\n\n")
for a in algorithms
  # print(" * Explicit ", uppercase(string(a)))
  # srand(seed); @time result = simulate(Simulation(m), t_final, a, o=Explicit(), itr=100_000)
  # computed = computed_mean(result)
  # @test_approx_eq_eps computed theoretical[end] 1e0

  print("   - Uniform ", uppercase(string(a)))
  srand(seed); @time result = simulate(m, tf=t_final, with=a, output=Uniform(), dt=0.1, itr=100_000)
  computed = aggregate(result, :Time, mean)
  print("     |observed - theoretical| = ", abs(computed[:X_mean][end] - theoretical[end]), "\n")
  println()
end
