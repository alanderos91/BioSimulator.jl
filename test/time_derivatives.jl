m = Simulation(kendall())
x = m.state
r = m.rxns
p = m.param

X = x[1]

α = p[:α]
μ = p[:μ]
ν = p[:ν]

# d/dt (x_1)
dxdt = [(1 * α * X) + (-1 * μ * X) + (1 * ν)]

# d/dt (r)
drdt = zeros(Float64, 3)
drdt[1] = BioSimulator.mass_action_deriv(r[1], x, p, 1) * dxdt[1]
drdt[2] = BioSimulator.mass_action_deriv(r[2], x, p, 1) * dxdt[1]
drdt[3] = BioSimulator.mass_action_deriv(r[3], x, p, 1) * dxdt[1]

# We have not computed propensities...
@test BioSimulator.compute_mean_derivatives!([0.0], r) == [0.0]
@test BioSimulator.compute_time_derivatives!([0.0, 0.0, 0.0], x, r, p, [0.0]) == [0.0, 0.0, 0.0]
@test BioSimulator.compute_time_derivatives!([1.0, 1.0, 1.0], x, r, p, [0.0]) == [0.0, 0.0, 0.0]

for rj in r; rj.propensity = BioSimulator.mass_action(rj, x, p); end

# We have computed propensities...
@test BioSimulator.compute_mean_derivatives!([0.0], r) == dxdt
@test BioSimulator.compute_time_derivatives!([0.0, 0.0, 0.0], x, r, p, dxdt) == drdt
@test BioSimulator.compute_time_derivatives!([1.0, 1.0, 1.0], x, r, p, dxdt) == drdt
