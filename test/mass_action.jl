m   = Simulation(test1())
X_t = m.state
r   = Dict{Symbol,BioSimulator.ReactionChannel}()
p   = m.param

x  = findfirst(m.sname, :X)
y  = findfirst(m.sname, :Y)
xy = findfirst(m.sname, :XY)

for reaction in m.rxns
  r[reaction.id] = reaction
end

print("   - Zero-Order Reaction:             ")
@test BioSimulator.mass_action(r[:zero], X_t, p)               == p[:k0]  # Propensity
@test BioSimulator.mass_action_deriv(r[:zero], X_t, p, x)      == 0.0     # Partial wrt X
@test BioSimulator.mass_action_deriv(r[:zero], X_t, p, y)      == 0.0     # Partial wrt Y
@test BioSimulator.mass_action_deriv(r[:zero], X_t, p, xy)     == 0.0     # Partial wrt XY
println("Passed")

print("   - First-Order Reaction:            ")
@test BioSimulator.mass_action(r[:first], X_t, p)              == p[:k1] * X_t[x]
@test BioSimulator.mass_action_deriv(r[:first], X_t, p, x)     == p[:k1]
@test BioSimulator.mass_action_deriv(r[:first], X_t, p, y)     == 0.0
@test BioSimulator.mass_action_deriv(r[:first], X_t, p, xy)    == 0.0
println("Passed")

print("   - Second-Order Reaction, Distinct: ")
@test BioSimulator.mass_action(r[:second_a], X_t, p)           == p[:k2a] * X_t[x] * X_t[y]
@test BioSimulator.mass_action_deriv(r[:second_a], X_t, p, x)  == p[:k2a] * X_t[y]
@test BioSimulator.mass_action_deriv(r[:second_a], X_t, p, y)  == p[:k2a] * X_t[x]
@test BioSimulator.mass_action_deriv(r[:second_a], X_t, p, xy) == 0.0
println("Passed")

print("   - Second-Order Reaction, Repeated: ")
@test BioSimulator.mass_action(r[:second_b], X_t, p)           == p[:k2b] * X_t[x] * (X_t[x] - 1)
@test BioSimulator.mass_action_deriv(r[:second_b], X_t, p, x)  == p[:k2b] * (2 * X_t[x] - 1)
@test BioSimulator.mass_action_deriv(r[:second_b], X_t, p, y)  == 0.0
@test BioSimulator.mass_action_deriv(r[:second_b], X_t, p, xy) == 0.0
println("Passed")
