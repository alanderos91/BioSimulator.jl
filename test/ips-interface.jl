import BioSimulator: IPSReactionStruct, DefaultIPSLaw

@testset "IPS" begin
  @testset "IPSReactionStruct" begin
    @testset "bad inputs" begin
      @testset "pairwise" begin
        # particle types cannot all be the same
        @test_throws ErrorException IPSReactionStruct(true, 1, 1, 1, 1, 2, 1.0, DefaultIPSLaw())

        # particle types must be positive
        @test_throws ErrorException IPSReactionStruct(true, 0, 1, 2, 3, 2, 1.0, DefaultIPSLaw())
        @test_throws ErrorException IPSReactionStruct(true, 1, 0, 2, 3, 2, 1.0, DefaultIPSLaw())
        @test_throws ErrorException IPSReactionStruct(true, 1, 2, 0, 3, 2, 1.0, DefaultIPSLaw())
        @test_throws ErrorException IPSReactionStruct(true, 1, 2, 3, 0, 2, 1.0, DefaultIPSLaw())

        # class must be a positive number
        @test_throws ErrorException IPSReactionStruct(true, 1, 2, 3, 4, 0, 1.0, DefaultIPSLaw())
        @test_throws ErrorException IPSReactionStruct(true, 1, 2, 3, 4, -1, 1.0, DefaultIPSLaw())

        # rate must be non-negative
        @test_throws ErrorException IPSReactionStruct(true, 1, 2, 3, 4, 5, -1.0, DefaultIPSLaw())
      end

      @testset "non-pairwise" begin
        # particle types cannot all be the same
        @test_throws ErrorException IPSReactionStruct(false, 1, 0, 1, 0, 2, 1.0, DefaultIPSLaw())

        # particle types must be positive
        @test_throws ErrorException IPSReactionStruct(false, 0, 0, 1, 0, 5, 1.0, DefaultIPSLaw())
        @test_throws ErrorException IPSReactionStruct(false, 1, 0, 0, 0, 5, 1.0, DefaultIPSLaw())

        # adjacent types must be zero
        @test_throws ErrorException IPSReactionStruct(false, 1, 0, 2, 1, 3, 1.0, DefaultIPSLaw())
        @test_throws ErrorException IPSReactionStruct(false, 1, 2, 2, 0, 3, 1.0, DefaultIPSLaw())
        @test_throws ErrorException IPSReactionStruct(false, 1, 2, 2, 1, 3, 1.0, DefaultIPSLaw())

        # class must be a positive number
        @test_throws ErrorException IPSReactionStruct(true, 1, 0, 2, 0, 0, 1.0, DefaultIPSLaw())
        @test_throws ErrorException IPSReactionStruct(true, 1, 0, 2, 0, -1, 1.0, DefaultIPSLaw())

        # rate must be non-negative
        @test_throws ErrorException IPSReactionStruct(true, 1, 0, 2, 0, 3, -1.0, DefaultIPSLaw())
      end
    end
    # need tests for checking that rates are computed correctly
  end

  @testset "macro interface" begin
    @testset "kinetic law parsing" begin
      # if no kinetic law is given, use DefaultIPSLaw
      intermediate = @def_reactions begin
        A + 0 --> A + A, α
      end α
      @test intermediate[1].klaw isa DefaultIPSLaw

      # use the specified kinetic law
      intermediate = @def_reactions begin
        A + 0 --> A + A, α, DefaultIPSLaw
      end α
      @test intermediate[1].klaw isa DefaultIPSLaw
    end

    @testset "bad inputs" begin
      # rate must be non-negative

      # empty sites must be the first type

      # assigned types must be positive

      # minimal number of classes
    end

      #= non-pairwise reactions
      A --> 0
      A --> ∅
      A --> B
      =#

      #= pairwise reactions
        A + 0 --> 0 + A
        A + ∅ --> ∅ + A
        A + ∅ --> A + A
        A + B --> A + A
        A + B --> A + C
        A + B --> C + D
      =#

    @testset "pure birth process" begin
      intermediate = @def_reactions begin
        A + 0 --> A + A, α
      end α

      params = [1.0]

      @testset "dimension = 1" begin
        model = @enumerate_with_sclass intermediate VonNeumann() 1 params

        rxn = model.reactions

        @test rxn[1] == IPSReactionStruct(true, 2, 1, 2, 2, 1, 1.0, DefaultIPSLaw())
        @test rxn[2] == IPSReactionStruct(true, 2, 1, 2, 2, 2, 2.0, DefaultIPSLaw())
      end

      @testset "dimension = 2" begin
        model = @enumerate_with_sclass intermediate VonNeumann() 2 params

        rxn = model.reactions

        @test rxn[1] == IPSReactionStruct(true, 2, 1, 2, 2, 1, 1.0, DefaultIPSLaw())
        @test rxn[2] == IPSReactionStruct(true, 2, 1, 2, 2, 2, 2.0, DefaultIPSLaw())
        @test rxn[3] == IPSReactionStruct(true, 2, 1, 2, 2, 3, 3.0, DefaultIPSLaw())
        @test rxn[4] == IPSReactionStruct(true, 2, 1, 2, 2, 4, 4.0, DefaultIPSLaw())
      end

      @testset "dimension = 3" begin
        model = @enumerate_with_sclass intermediate VonNeumann() 3 params

        rxn = model.reactions

        @test rxn[1] == IPSReactionStruct(true, 2, 1, 2, 2, 1, 1.0, DefaultIPSLaw())
        @test rxn[2] == IPSReactionStruct(true, 2, 1, 2, 2, 2, 2.0, DefaultIPSLaw())
        @test rxn[3] == IPSReactionStruct(true, 2, 1, 2, 2, 3, 3.0, DefaultIPSLaw())
        @test rxn[4] == IPSReactionStruct(true, 2, 1, 2, 2, 4, 4.0, DefaultIPSLaw())
        @test rxn[5] == IPSReactionStruct(true, 2, 1, 2, 2, 5, 5.0, DefaultIPSLaw())
        @test rxn[6] == IPSReactionStruct(true, 2, 1, 2, 2, 6, 6.0, DefaultIPSLaw())
      end
    end

    @testset "pure death process" begin
      intermediate = @def_reactions begin
        X --> 0, μ
      end μ

      params = [1.0]

      @testset "dimension = 1" begin
        model = @enumerate_with_sclass intermediate VonNeumann() 1 params

        rxn = model.reactions

        @test rxn[1] == IPSReactionStruct(false, 2, 0, 1, 0, 1, 1.0, DefaultIPSLaw())
      end

      @testset "dimension = 2" begin
        model = @enumerate_with_sclass intermediate VonNeumann() 2 params

        rxn = model.reactions

        @test rxn[1] == IPSReactionStruct(false, 2, 0, 1, 0, 1, 1.0, DefaultIPSLaw())
      end

      @testset "dimension = 3" begin
        model = @enumerate_with_sclass intermediate VonNeumann() 3 params

        rxn = model.reactions

        @test rxn[1] == IPSReactionStruct(false, 2, 0, 1, 0, 1, 1.0, DefaultIPSLaw())
      end
    end

    @testset "pure migration process" begin
      intermediate = @def_reactions begin
        X + 0 --> 0 + X, ν
      end ν

      params = [1.0]

      @testset "dimension = 1" begin
        model = @enumerate_with_sclass intermediate VonNeumann() 1 params

        rxn = model.reactions

        @test rxn[1] == IPSReactionStruct(true, 2, 1, 1, 2, 1, 1.0, DefaultIPSLaw())
        @test rxn[2] == IPSReactionStruct(true, 2, 1, 1, 2, 2, 2.0, DefaultIPSLaw())
      end

      @testset "dimension = 2" begin
        model = @enumerate_with_sclass intermediate VonNeumann() 2 params

        rxn = model.reactions

        @test rxn[1] == IPSReactionStruct(true, 2, 1, 1, 2, 1, 1.0, DefaultIPSLaw())
        @test rxn[2] == IPSReactionStruct(true, 2, 1, 1, 2, 2, 2.0, DefaultIPSLaw())
        @test rxn[3] == IPSReactionStruct(true, 2, 1, 1, 2, 3, 3.0, DefaultIPSLaw())
        @test rxn[4] == IPSReactionStruct(true, 2, 1, 1, 2, 4, 4.0, DefaultIPSLaw())
      end

      @testset "dimension = 3" begin
        model = @enumerate_with_sclass intermediate VonNeumann() 3 params

        rxn = model.reactions

        @test rxn[1] == IPSReactionStruct(true, 2, 1, 1, 2, 1, 1.0, DefaultIPSLaw())
        @test rxn[2] == IPSReactionStruct(true, 2, 1, 1, 2, 2, 2.0, DefaultIPSLaw())
        @test rxn[3] == IPSReactionStruct(true, 2, 1, 1, 2, 3, 3.0, DefaultIPSLaw())
        @test rxn[4] == IPSReactionStruct(true, 2, 1, 1, 2, 4, 4.0, DefaultIPSLaw())
        @test rxn[5] == IPSReactionStruct(true, 2, 1, 1, 2, 5, 5.0, DefaultIPSLaw())
        @test rxn[6] == IPSReactionStruct(true, 2, 1, 1, 2, 6, 6.0, DefaultIPSLaw())
      end
    end
  end
end
