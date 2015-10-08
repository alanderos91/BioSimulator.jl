# Kendall's Process
x = [Species("Particle", 5, true)]
r = [Reaction("Birth", "alpha", [1], [2]),
     Reaction("Death", "mu", [1], [0]),
     Reaction("Immigration", "nu", [0], [1])]
p = Dict{ASCIIString, Float64}("alpha" => 2.0,
                               "mu" => 1.0,
                               "nu" => 0.5)
kendall = Network("Kendall's Process", x, r, p);

# Two-Type Branching Process
x = [Species("Male", 5, true),
     Species("Female", 5, true)]
r = [Reaction("Female Birth", "lambdaP", [1,0], [2,0]),
     Reaction("Male Birth", "lambdaQ", [1,0], [1,1]),
     Reaction("Female Death", "mu", [1,0], [0,0]),
     Reaction("Male Death", "nu", [0,1], [0,0])]
p = Dict{ASCIIString, Float64}("lambdaP" => 4.0 * 0.45,
                               "lambdaQ" => 4.0 * 0.55,
                               "mu" => 1.0,
                               "nu" => 1.0)
twotype = Network("Two-Type Branching Process", x, r, p);

# Michaelis-Menten Enzyme Kinetics
x = [Species("S1", 300, true),
     Species("S2", 120, true),
     Species("S3", 0, true),
     Species("S4", 0, true)]
r = [Reaction("Dimerization", "a1", [1,1,0,0], [0,0,1,0]),
     Reaction("Dissociation", "a2", [0,0,1,0], [1,1,0,0]),
     Reaction("Conversion", "a3", [0,0,1,0], [0,1,0,1])]
p = Dict{ASCIIString, Float64}("a1" => 0.00166,
                               "a2" => 0.001,
                               "a3" => 0.1)
michment = Network("Michaelis Menten Enzyme Kinetics", x, r, p);

# Simulation Parameters
tf = 4.0
dt = 4.0
seed = 5357

# SSA Convergence Tests
srand(seed)
result = ssa(Simulation(kendall), tf, itr=10^5)
computed = mean(result)
for i in eachindex(theoretical)
  @test_approx_eq_eps mean(result)[i] theoretical[i] 1.0
end

srand(seed)
result = ssa(Simulation(twotype), tf, itr=10^5)
computed = mean(result)
for i in eachindex(theoretical)
  @test_approx_eq_eps mean(result)[i] theoretical[i] 1.0
end

# dSSA Convergence Tests
srand(seed)
result = dssa(Simulation(kendall), tf, dt=dt itr=10^5)
computed = mean(result)
for i in eachindex(theoretical)
  @test_approx_eq_eps computed[i] theoretical[i] 1.0
end

srand(seed)
result = dssa(Simulation(twotype), tf, dt=dt itr=10^5)
computed = mean(result)
for i in eachindex(theoretical)
  @test_approx_eq_eps computed[i] theoretical[i] 1.0
end
