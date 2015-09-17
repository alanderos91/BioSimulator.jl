# Containers
X_t = Species[]
r = Dict{Symbol, Reaction}()
param = Dict{ASCIIString, Float64}()

# Species
push!(X_t, Species("X", 100, true)) # 1
push!(X_t, Species("Y", 500, true)) # 2
push!(X_t, Species("XY",  1, true)) # 3

# Zero-Order Reaction: 0 ---> X
r[:zero] = Reaction("Immigration", "k0", 0.0, [0,0,0], [1,0,0])
param["k0"] = 0.1

# First-Order Reaction: X ---> 0
r[:first] = Reaction("Death", "k1", 0.0, [1,0,0],[0,0,0])
param["k1"] = 0.01

# Second-Order Reaction: X + Y ---> XY
r[:second_a] = Reaction("Dimerization_XY", "k2a", 0.0, [1,1,0], [0,0,1])
param["k2a"] = 0.25

# Second-Order Reaction: X + X ---> Y
r[:second_b] = Reaction("Dimerization_Y", "k2b", 0.0, [2,0,0], [0,1,0])
param["k2b"] = 0.5

@test BioSimulator.mass_action(r[:zero], X_t, param) == param["k0"]  # Propensity
@test BioSimulator.mass_action_deriv(r[:zero], X_t, param, 1) == 0.0 # Partial wrt X
@test BioSimulator.mass_action_deriv(r[:zero], X_t, param, 2) == 0.0 # Partial wrt Y
@test BioSimulator.mass_action_deriv(r[:zero], X_t, param, 3) == 0.0 # Partial wrt XY

@test BioSimulator.mass_action(r[:first], X_t, param) == param["k1"] * X_t[1].pop
@test BioSimulator.mass_action_deriv(r[:first], X_t, param, 1) == param["k1"]
@test BioSimulator.mass_action_deriv(r[:first], X_t, param, 2) == 0.0
@test BioSimulator.mass_action_deriv(r[:first], X_t, param, 3) == 0.0

@test BioSimulator.mass_action(r[:second_a], X_t, param) == param["k2a"] * X_t[1].pop * X_t[2].pop
@test BioSimulator.mass_action_deriv(r[:second_a], X_t, param, 1) == param["k2a"] * X_t[2].pop
@test BioSimulator.mass_action_deriv(r[:second_a], X_t, param, 2) == param["k2a"] * X_t[1].pop
@test BioSimulator.mass_action_deriv(r[:second_a], X_t, param, 3) == 0.0

@test BioSimulator.mass_action(r[:second_b], X_t, param) == param["k2b"] * X_t[1].pop * (X_t[1].pop - 1)
@test BioSimulator.mass_action_deriv(r[:second_b], X_t, param, 1) == param["k2b"] * (2 * X_t[1].pop - 1)
@test BioSimulator.mass_action_deriv(r[:second_b], X_t, param, 2) == 0.0
@test BioSimulator.mass_action_deriv(r[:second_b], X_t, param, 3) == 0.0
