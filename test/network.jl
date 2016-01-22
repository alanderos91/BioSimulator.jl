m = Network("Test")

# Should not add Reaction to model without Species

@test_throws ErrorException m <= Reaction(:Bad, :k0, r=(:X => 1))
@test_throws ErrorException m <= Reaction(:Bad, :k0, p=(:X => 1))
@test length(m.reactions) == 0

# Should not add Reaction with undefined Species
m <= Species(:X, 100)
m <= Species(:Y, 100)

@test_throws ErrorException m <= Reaction(:Bad, :k0, r=(:X => 1), p=(:Z => 1))
@test_throws ErrorException m <= Reaction(:Bad, :k0, p=(:Z => 1), r=(:Y => 1))
@test length(m.reactions) == 0
