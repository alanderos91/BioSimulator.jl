using BioSimulator
using Base.Test
using DataFrames

function kendall_mean(i,t,alpha,mu,nu)
   x = exp((alpha-mu)*t)
  return i * x + nu/(alpha-mu)*(x-1)
end

function computed_mean(sj)
  mval = 0.0
  for r in sj.results
    mval = mval + r["X"].states[end].pop
  end
  return mval / length(sj.results)
end

# Values
x = 5
α = 2.0
μ = 1.0
ν = 0.5

# Model Definition
network = Network("Kendall's Process")

add_species!(network, "X", x)

add_reaction!(network, "Birth", "alpha")
add_reactant!(network, "Birth", "X", 1)
add_product!(network, "Birth", "X", 2)
set_parameter!(network, "alpha", α)

add_reaction!(network, "Death", "mu")
add_reactant!(network, "Death", "X")
set_parameter!(network, "mu", μ)

add_reaction!(network, "Immigration", "nu")
add_product!(network, "Immigration", "X")
set_parameter!(network, "nu", ν)

t_final = 4.0
Δt = 0.1
seed = 5357

# Compute mean for comparison
points = round(Int, t_final/Δt) + 1
t = linspace(0.0,t_final, points)
theoretical = kendall_mean(x,t,α,μ,ν)

algorithms = [:ssa, :odm, :nrm, :sal, :frm]

# Run SSA and SAL once to compile
print("Precompiling..."); @time begin
  for a in algorithms
    simulate(network, tf=t_final, with=a, output=Explicit(), itr=1)
    simulate(network, tf=t_final, with=a, output=Uniform(), dt=0.1, itr=1)
  end
end

println("Running tests...")
for a in algorithms
  # print(" * Explicit ", uppercase(string(a)))
  # srand(seed); @time result = simulate(Simulation(network), t_final, a, o=Explicit(), itr=100_000)
  # computed = computed_mean(result)
  # @test_approx_eq_eps computed theoretical[end] 1e0

  print(" *  Uniform ", uppercase(string(a)))
  srand(seed); @time result = simulate(network, tf=t_final, with=a, output=Uniform(), dt=0.1, itr=100_000)
  computed = aggregate(result, :Time, mean)
  @test_approx_eq_eps computed[:X_mean][end] theoretical[end] 1e0
end
