network = test1()
Xt, species_id, id2ind = BioSimulator.make_species_arr(network.species)
rv, reaction_id        = BioSimulator.make_reaction_arr(network.reactions, id2ind)

m  = BioSimulator.Model(network.id, Xt, rv, network.parameters, deepcopy(Xt))
p  = BioSimulator.parameters(m)
r  = Dict{Symbol,BioSimulator.ReactionChannel}()
x  = findfirst(species_id, :X)
y  = findfirst(species_id, :Y)
xy = findfirst(species_id, :XY)

for reaction in rv.reactions
  r[reaction.id] = reaction
end

print("   - Zero-Order Reaction:             ")
@test BioSimulator.mass_action(r[:zero], Xt, p)               == p[:k0]  # Propensity
@test BioSimulator.mass_action_deriv(r[:zero], Xt, p, x)      == 0.0     # Partial wrt X
@test BioSimulator.mass_action_deriv(r[:zero], Xt, p, y)      == 0.0     # Partial wrt Y
@test BioSimulator.mass_action_deriv(r[:zero], Xt, p, xy)     == 0.0     # Partial wrt XY
println("Passed")

print("   - First-Order Reaction:            ")
@test BioSimulator.mass_action(r[:first], Xt, p)              == p[:k1] * Xt[x]
@test BioSimulator.mass_action_deriv(r[:first], Xt, p, x)     == p[:k1]
@test BioSimulator.mass_action_deriv(r[:first], Xt, p, y)     == 0.0
@test BioSimulator.mass_action_deriv(r[:first], Xt, p, xy)    == 0.0
println("Passed")

print("   - Second-Order Reaction, Distinct: ")
@test BioSimulator.mass_action(r[:second_a], Xt, p)           == p[:k2a] * Xt[x] * Xt[y]
@test BioSimulator.mass_action_deriv(r[:second_a], Xt, p, x)  == p[:k2a] * Xt[y]
@test BioSimulator.mass_action_deriv(r[:second_a], Xt, p, y)  == p[:k2a] * Xt[x]
@test BioSimulator.mass_action_deriv(r[:second_a], Xt, p, xy) == 0.0
println("Passed")

print("   - Second-Order Reaction, Repeated: ")
@test BioSimulator.mass_action(r[:second_b], Xt, p)           == p[:k2b] * Xt[x] * (Xt[x] - 1)
@test BioSimulator.mass_action_deriv(r[:second_b], Xt, p, x)  == p[:k2b] * (2 * Xt[x] - 1)
@test BioSimulator.mass_action_deriv(r[:second_b], Xt, p, y)  == 0.0
@test BioSimulator.mass_action_deriv(r[:second_b], Xt, p, xy) == 0.0
println("Passed")
