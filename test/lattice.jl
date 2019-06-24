import NewBioSimulator: get_ptype, get_neighbor_class, initialize_datastructs!, coordinates, label

@testset "Lattice" begin
  @testset "initialization" begin
    ##### reaction definition #####
    intermediate = @def_reactions begin
      A + B --> A + A, k1 # predation / predator reproduction
      B + 0 --> 0 + B, k2 # prey movement
      A + 0 --> 0 + A, k3 # predator movement
      B + 0 --> B + B, k4 # prey reproduction
      A --> 0, k5         # predator death
      B --> 0, k6         # prey death
    end k1 k2 k2 k3 k4 k5 k6

    ##### assign parameter values #####
    params = [2.0, 1.0, 1.0, 2.0, 1.0, 1.0]

    ##### generate the full reaction list #####
    model = @enumerate_with_sclass intermediate VonNeumann() 2 params

    ##### initial configuration #####
    initial_sites = [(0,1), (1,0), (1,1), (1,2), (1,3), (2,0), (2,1), (2,2), (3,0), (3,1)]
    initial_sites = [initial_sites[j][i] for i in 1:2, j in 1:length(initial_sites)]
    ptypes        = [1, 1, 1, 2, 2, 1, 2, 2, 2, 2]

    lattice = Lattice(initial_sites, ptypes)

    ##### correct answers #####
    expected_opensites = [
    (0, 0), (1, -1), (2, -1), (3, -1), (4, 0),
    (4, 1), (3, 2), (2, 3), (1, 4), (0, 3),
    (0, 2), (-1, 1)
    ]

    expected_compositions = [
    [4, 0, 0], [3, 1, 0], [3, 0, 1], [2, 2, 0], [2, 1, 1],
    [2, 0, 2], [1, 3, 0], [1, 2, 1], [1, 1, 2], [1, 0, 3],
    [0, 4, 0], [0, 3, 1], [0, 2, 2], [0, 1, 3], [0, 0, 4]
    ]

    expected_matrix = zeros(Int, 15, 3)

    # open site counts
    expected_matrix[2, 1] = 3
    expected_matrix[3, 1] = 5
    expected_matrix[4, 1] = 1
    expected_matrix[5, 1] = 1
    expected_matrix[6, 1] = 2
    # total = 12

    # predator counts
    expected_matrix[2,  2] = 1
    expected_matrix[4,  2] = 1
    expected_matrix[9,  2] = 1
    expected_matrix[13, 2] = 1
    # total = 4

    # prey counts
    expected_matrix[3,  3] = 1
    expected_matrix[5,  3] = 1
    expected_matrix[6,  3] = 2
    expected_matrix[9,  3] = 1
    expected_matrix[13, 3] = 1
    # total = 6

    expected_neighborclass = Dict(
    # particles
    (0, 1) => 2,
    (1, 0) => 4,
    (1, 1) => 13,
    (1, 2) => 9,
    (1, 3) => 3,
    (2, 0) => 9,
    (2, 1) => 13,
    (2, 2) => 6,
    (3, 0) => 5,
    (3, 1) => 6,
    # open sites
    (0, 0) => 4,
    (1, -1) => 2,
    (2, -1) => 2,
    (3, -1) => 3,
    (4, 0) => 3,
    (4, 1) => 3,
    (3, 2) => 6,
    (2, 3) => 6,
    (1, 4) => 3,
    (0, 3) => 3,
    (0, 2) => 5,
    (-1, 1) => 2
    )
    expected_sampleclass = Dict(
    # particles only
    (0, 1) => [11, 13],
    (1, 0) => [10, 13],
    (1, 1) => [2, 13],
    (1, 2) => [5, 14],
    (1, 3) => [7, 14],
    (2, 0) => [9, 2, 13],
    (2, 1) => [14],
    (2, 2) => [6, 14],
    (3, 0) => [6, 14],
    (3, 1) => [6, 14]
    )

    # perform the initialization step
    initialize_datastructs!(lattice, model)

    @testset "open site tracking" begin
      test_opensites = filter(x -> isequal(get_ptype(x), 1), lattice.site)

      # do we generate the correct number of open sites?
      @test length(test_opensites) == length(expected_opensites)

      # do we spawn the correct open sites?
      test_coordinates = [tuple(coordinates(x)...) for x in lattice.site]

      for coord in expected_opensites
        @test coord in test_coordinates
      end
    end

    @testset "compositions" begin
      test_compositions = model.enumeration.composition

      # do we generate the coorect number of compositions?
      @test length(test_compositions) == length(expected_compositions)

      # do we generate the correct compositions?
      for composition in expected_compositions
        @test composition in test_compositions
      end

      # do we iterate through the compositions in a specific order (open sites -> particle types...)?
      for k in eachindex(test_compositions)
        composition = test_compositions[k]
        @test composition == expected_compositions[k]
      end
    end

    @testset "neighborhood classes" begin
      # is each site assigned to the correct neighborhood class?
      for x in lattice.site
        coord = tuple(coordinates(x)...)
        @test get_neighbor_class(x) == expected_neighborclass[coord]
      end
    end

    @testset "sample classes" begin
      test_sampleclass = model.enumeration.class

      # is each site assigned to the correct sample classes?
      for x in lattice.site
        if get_ptype(x) > 1
          coord = tuple(coordinates(x)...)
          for s in expected_sampleclass[coord]
            @test label(x) in test_sampleclass[s]
          end
        end
      end
    end

    @testset "count matrix" begin
      test_matrix = zero(expected_matrix)

      for x in lattice.site
        l = get_ptype(x)
        k = get_neighbor_class(x)

        test_matrix[k,l] += 1
      end

      # is the count matrix correct?
      for i in eachindex(test_matrix)
        @test test_matrix[i] == expected_matrix[i]
      end
    end
  end
end
