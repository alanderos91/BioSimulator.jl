network = kendall()
Xt, species_id, id2ind = BioSimulator.make_species_arr(network.species)
rv, reaction_id        = BioSimulator.make_reaction_arr(network.reactions, id2ind)

m = BioSimulator.Model(network.id, Xt, rv, network.parameters, deepcopy(Xt))
rv = BioSimulator.reactions(m)
p = BioSimulator.parameters(m)

X = Xt[1]

α = p[:α]
μ = p[:μ]
ν = p[:ν]

# d/dt (x_1)
dxdt = [(1 * α * X) + (-1 * μ * X) + (1 * ν)]

# d/dt (r)
drdt = zeros(Float64, 3)
drdt[1] = BioSimulator.mass_action_deriv(BioSimulator.reaction(rv, 1), Xt, p, 1) * dxdt[1]
drdt[2] = BioSimulator.mass_action_deriv(BioSimulator.reaction(rv, 2), Xt, p, 1) * dxdt[1]
drdt[3] = BioSimulator.mass_action_deriv(BioSimulator.reaction(rv, 3), Xt, p, 1) * dxdt[1]

# We have not computed propensities...
@test BioSimulator.mean_derivatives!([0.0], m) == [0.0]
@test BioSimulator.time_derivatives!([0.0, 0.0, 0.0], m, [0.0, 0.0, 0.0]) == [0.0, 0.0, 0.0]
@test BioSimulator.time_derivatives!([1.0, 1.0, 1.0], m, [0.0, 0.0, 0.0]) == [0.0, 0.0, 0.0]

BioSimulator.compute_propensities!(m)

# We have computed propensities...
@test BioSimulator.mean_derivatives!([0.0], m) == dxdt
@test BioSimulator.time_derivatives!([0.0, 0.0, 0.0], m, dxdt) == drdt
@test BioSimulator.time_derivatives!([1.0, 1.0, 1.0], m, dxdt) == drdt
