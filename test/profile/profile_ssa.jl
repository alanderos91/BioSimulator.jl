using BioSimulator

x = [Species("Particle", 5, true)]
r = [Reaction("Birth", "alpha", [1], [2]),
     Reaction("Death", "mu", [1], [0]),
     Reaction("Immigration", "nu", [0], [1])]
p = Dict{ASCIIString, Float64}("alpha" => 2.0,
                               "mu" => 1.0,
                               "nu" => 0.5)
kendall = Network("Kendall's Process", x, r, p);

# Compile
@time ssa(Simulation(kendall), 4.0)
Profile.clear_malloc_data()

# Test
@time ssa(Simulation(kendall), 4.0)
