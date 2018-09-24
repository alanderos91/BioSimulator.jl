@testset "Network" begin
  m = Network("Test")

  # Adding a Reaction with an unknown Species should throw an error
  @test_throws ErrorException m <= Reaction("bad",      1.0, "X --> 0")
  @test_throws ErrorException m <= Reaction("also bad", 1.0, "0 --> X")
  @test number_reactions(m) == 0

  # Adding a Species with negative counts should throw an error
  @test_throws ErrorException m <= Species("X", -1)

  # Adding a Species with non-negative counts is okay
  m <= Species("X", 0)
  m <= Species("Y", 100)

  @test number_species(m) == 2

  # Adding a Reaction with an unknown Species should throw an error
  @test_throws ErrorException m <= Reaction("should not work",  1.0, "X --> Z")
  @test_throws ErrorException m <= Reaction("should also fail", 1.0, "Z --> Y")
  @test number_reactions(m) == 0

  # Adding the same Species twice should change its initial value
  m <= Species("X", 50)

  @test number_species(m) == 2
  @test get_species(m, "X").population == 50
end