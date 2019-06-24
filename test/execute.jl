import NewBioSimulator: execute_jump!, get_site, initialize_datastructs!, get_neighbor_class, get_ptype, coordinates

@testset "Execute" begin
  @testset "Predator-Prey Model" begin
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

    initialize_datastructs!(lattice, model)
    enum = model.enumeration

    # Test Reaction:
    execute_jump!(lattice, model, 9)

    # aggregate coordinates
    test_coord = [tuple(coordinates(x)...) for x in lattice.site]

    # check new open sites
    @test (2, -1) in test_coord

    # check type changes
    @test get_ptype(get_site(lattice, (2,  0))) == 1
    @test get_ptype(get_site(lattice, (2, -1))) == 2

    # check neighbor class changes
    # (2, 1): 13 -> 9
    # (1, 0): 4  -> 2
    # (3, 0): 5  -> 3
    # (2, 0): 9  -> 13
    # (2, -1): 2 -> 1
    points_to_check = [(2,1), (1,0), (3,0), (2,0), (2,-1)]
    list = [get_site(lattice, x) for x in points_to_check]

    @test get_neighbor_class(list[1]) == 9
    @test get_neighbor_class(list[2]) == 2
    @test get_neighbor_class(list[3]) == 3
    @test get_neighbor_class(list[4]) == 13
    @test get_neighbor_class(list[5]) == 1

    # check sample class changes
    # (2, 1): [14] -> [14, 5]
    # (1, 0): [10, 13] -> [13, 11]
    # (3, 0): [6, 14] -> [14, 7]
    # (2, 0): [2,9,13] -> []
    # (2, -1): [] -> [12, 13]
    particle_sample_classes = [
      (2, 1)  => [5, 14],
      (1, 0)  => [11, 13],
      (3, 0)  => [7, 14],
      (2, 0)  => Int[],
      (2, -1) => [12, 13]
      ]

    for (coord, expected_sample_classes) in particle_sample_classes
      x = get_site(lattice, coord)

      l = get_ptype(x)
      k = get_neighbor_class(x)

      test_sample_classes = sort(get(enum.pair_to_classes, (l,k), Int[]))

      @test test_sample_classes == expected_sample_classes
    end

    # check invariants
    @test issorted(lattice.site)
    @test issorted(lattice.coord_order, by = coordinates)
  end
end
