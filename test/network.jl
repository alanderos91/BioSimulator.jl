m = Network("Test")

# Should not add Reaction to model without Species

@test_throws ErrorException m <= Reaction(:Bad, :k0, :(X --> 0))
@test_throws ErrorException m <= Reaction(:Bad, :k0, :(0 --> X))
@test length(m.reactions) == 0

# Should not add Reaction with undefined Species
m <= Species(:X, 100)
m <= Species(:Y, 100)

@test_throws ErrorException m <= Reaction(:Bad, :k0, :(X --> Z))
@test_throws ErrorException m <= Reaction(:Bad, :k0, :(Z --> Y))
@test length(m.reactions) == 0
