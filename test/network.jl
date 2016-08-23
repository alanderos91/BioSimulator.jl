m = Network("Test")

# Should not add Reaction to model without Species
@test_throws ErrorException m <= Reaction(:Bad, 1.0, :(X --> 0))
@test_throws ErrorException m <= Reaction(:Bad, 1.0, :(0 --> X))
@test n_reactions(m) == 0

# Should not add Reaction with undefined Species
m <= Species(:X, 100)
m <= Species(:Y, 100)

@test_throws ErrorException m <= Reaction(:Bad, 1.0, :(X --> Z))
@test_throws ErrorException m <= Reaction(:Bad, 1.0, :(Z --> Y))
@test n_reactions(m) == 0

# Adding the same Species should change its initial value
m <= Species(:X, 50)
@test n_species(m) == 2
@test species_list(m)[:X].population == 50

# # Remove a Species
# m >= (:species, UTF8String("Y"))
# @test n_species(m) == 1
