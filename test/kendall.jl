using BioSimulator
using Base.Test

function kendall_mean(i,t,alpha,mu,nu)
   x = exp((alpha-mu)*t)
  return i * x + nu/(alpha-mu)*(x-1)
end

function computed_mean(sr,itr)
  mval = 0.0
  for i = 1:itr
    mval = mval + sr.results[i]["Particle"].states[end].value
  end
  return mval / itr
end

x = [Species("Particle", 5, true)]

r1 = Reaction("Birth", "alpha", [1], [2])
r2 = Reaction("Death", "mu", [1], [0])
r3 = Reaction("Immigration", "nu", [0], [1])
r = [r1, r2, r3]

p = Dict{ASCIIString, Float64}("alpha" => 2.0, "mu" => 1.0, "nu" => 0.5)
network = Network("Kendall's Process", x, r, p);

t_final = 4.0
Δt = 0.1

# Compute mean for comparison
points = round(Int, t_final/Δt) + 1
t = linspace(0.0,t_final, points)
m = kendall_mean(x[1].pop,t,p["alpha"],p["mu"],p["nu"])

# Run SSA and SAL once to compile
srand(5137)
ssa1 = ssa(Simulation(network), t_final)
srand(5137)
sal1 = sal(Simulation(network), t_final)

# Run model with SSA for 10^5 iterations
srand(5137)
@time ssa2 = ssa(Simulation(network), t_final, tracing=true, output=:fixed, stepsize=Δt, itr=10^5)
@test_approx_eq_eps computed_mean(ssa2,10^5) m[end] 2.0

# Run model with SAL for 10^5 iterations
srand(5137)
@time sal2 = sal(Simulation(network), t_final, tracing=true, output=:fixed, stepsize=Δt, itr=10^5)
@test_approx_eq_eps computed_mean(sal2,10^5) m[end] 2.0
