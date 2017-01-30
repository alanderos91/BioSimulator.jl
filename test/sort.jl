import BioSimulator: stoichiometry, coefficients, scaled_rates, propensities, dependencies
m = autoreg()

c = n_species(m)
d = n_reactions(m)

species   = species_list(m)
reactions = reaction_list(m)

X0, id, id2ind = BioSimulator.make_species_vector(species)

# dense
r1 = BioSimulator.DenseReactionSystem(reactions, id2ind, c, d)

# test 1
r2 = r1

@test stoichiometry(r1) == stoichiometry(r2)
@test coefficients(r1) == coefficients(r2)
@test scaled_rates(r1) == scaled_rates(r2)
@test dependencies(r1) == dependencies(r2)

# test 2
r2 = deepcopy(r1)

reaction_events = [1, 2, 3, 4, 5, 6, 7, 8]
ix = sortperm(reaction_events)
BioSimulator._sort!(r2, ix)

@test stoichiometry(r1) == stoichiometry(r2)
@test coefficients(r1) == coefficients(r2)
@test scaled_rates(r1) == scaled_rates(r2)
@test dependencies(r1) == dependencies(r2)

# test 3
r2 = deepcopy(r1)
reaction_events = [8, 7, 6, 5, 4, 3, 2, 1]
ix = sortperm(reaction_events)
BioSimulator._sort!(r2, ix)

@test stoichiometry(r1) != stoichiometry(r2)
@test coefficients(r1) != coefficients(r2)
@test scaled_rates(r1) != scaled_rates(r2)
@test dependencies(r1) != dependencies(r2)

# # sparse
# r = SparseReactionSystem(reactions, id2ind, c, d)
