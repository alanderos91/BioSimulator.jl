x = [Species("Particle", 5, true)]

r1 = Reaction("Birth", "alpha", [1], [2])
r2 = Reaction("Death", "mu", [1], [0])
r3 = Reaction("Immigration", "nu", [0], [1])
r = [r1, r2, r3]

p = Dict{ASCIIString, Float64}("alpha" => 2.0, "mu" => 1.0, "nu" => 0.5)
network = Network("Kendall's Process", x, r, p);

# d/dt (x_1)
dxdt_1 = (1 * 2.0 * 5) + (-1 * 1.0 * 5) + (1 * 0.5)

# d/dt (r)
drdt = zeros(Float64, 3)
drdt[1] = BioSimulator.mass_action_deriv(r[1], x, p, 1) * dxdt_1
drdt[2] = BioSimulator.mass_action_deriv(r[2], x, p, 1) * dxdt_1
drdt[3] = BioSimulator.mass_action_deriv(r[3], x, p, 1) * dxdt_1

# We have not computed propensities...
@test BioSimulator.mean_derivative(r, 1) == 0.0
@test BioSimulator.compute_time_derivatives!([0.0, 0.0, 0.0], x, r, p) == [0.0, 0.0, 0.0]
@test BioSimulator.compute_time_derivatives!([1.0, 1.0, 1.0], x, r, p) == [0.0, 0.0, 0.0]

for rj in r; rj.propensity = BioSimulator.mass_action(rj, x, p); end

# We have computed propensities...
@test BioSimulator.mean_derivative(r, 1) == dxdt_1
@test BioSimulator.compute_time_derivatives!([0.0, 0.0, 0.0], x, r, p) == drdt
@test BioSimulator.compute_time_derivatives!([1.0, 1.0, 1.0], x, r, p) == drdt
