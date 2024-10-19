import BioSimulator: MassAction
import BioSimulator: ReactionStruct, execute_jump!, rate, rate_derivative
import BioSimulator: ReactionSystem

# define notion of equality for ReactionStruct
import Base: ==

function ==(r1::ReactionStruct{T}, r2::ReactionStruct{S}) where {T,S}
  isequal = true
  isequal = isequal && T == S
  isequal = isequal && length(r1.reactants) == length(r2.reactants)
  isequal = isequal && all(x -> x in r2.reactants, r1.reactants)
  isequal = isequal && length(r1.net_change) == length(r2.net_change)
  isequal = isequal && all(x -> x in r2.net_change, r2.net_change)

  return isequal
end

@testset "ReactionSystem" begin
  x = rand(0:1000, 5)
  params = 10 * rand(4)

  # Species | A, B, C, D E
  #   index | 1  2  3  4 5

  # Order 0
  order0 = [
    (Tuple{Int,Int}[], [(1, 1)], 1), # 0 --> A
    (Tuple{Int,Int}[], [(2, 4)], 1)  # 0 --> 4 B
  ]

  # Order 1
  order1 = [
    ([(1, 1)], [(1, -1)],                 2), # A --> 0
    ([(1, 1)], [(1, -1), (2, 1)],         2), # A --> B
    ([(1, 1)], [(1, -1), (3, 2), (5, 2)], 2)  # A --> 2 C + 2 E
  ]

  # Order 2A
  order2a = [
    ([(1, 1), (2, 1)], [(1, -1), (2, -1), (3, 1)], 3) # A + B --> C
  ]

  # Order 2B
  order2b = [
    ([(2, 2)], [(2, -2), (1, 1)], 4) # 2 B --> A
  ]

  # Order 3
  # A + B + C --> D
  order3a = [
    ([(1, 1), (2, 1), (3, 1)], [(1, -1), (2, -1), (3, -1), (4, 1)], 1)
  ]

  # 3 * A + B + 2 * C --> D
  order3b = [
    ([(1, 3), (2, 1), (3, 2)], [(1, -3), (2, -1), (3, -2), (4, 1)], 1)
  ]

  @testset "ReactionStruct" begin
    rxn0  = [ReactionStruct(MassAction(), r, v, i) for (r, v, i) in order0]
    rxn1  = [ReactionStruct(MassAction(), r, v, i) for (r, v, i) in order1]
    rxn2a = [ReactionStruct(MassAction(), r, v, i) for (r, v, i) in order2a]
    rxn2b = [ReactionStruct(MassAction(), r, v, i) for (r, v, i) in order2b]
    rxn3a = [ReactionStruct(MassAction(), r, v, i) for (r, v, i) in order3a]
    rxn3b = [ReactionStruct(MassAction(), r, v, i) for (r, v, i) in order3b]

    @testset "fire_reaction!" begin
      result   = copy(x)
      expected = x + [1, 0, 0, 0, 0]
      execute_jump!(result, rxn0[1])
      @test result == expected

      result   = copy(x)
      expected = x + [0, 4, 0, 0, 0]
      execute_jump!(result, rxn0[2])
      @test result == expected

      result   = copy(x)
      expected = x + [-1, 0, 0, 0, 0]
      execute_jump!(result, rxn1[1])
      @test result == expected

      result   = copy(x)
      expected = x + [-1, 1, 0, 0, 0]
      execute_jump!(result, rxn1[2])
      @test result == expected

      result   = copy(x)
      expected = x + [-1, 0, 2, 0, 2]
      execute_jump!(result, rxn1[3])
      @test result == expected

      result   = copy(x)
      expected = x + [-1, -1, 1, 0, 0]
      execute_jump!(result, rxn2a[1])
      @test result == expected

      result   = copy(x)
      expected = x + [1, -2, 0, 0, 0]
      execute_jump!(result, rxn2b[1])
      @test result == expected

      result   = copy(x)
      expected = x + [-1, -1, -1, 1, 0]
      execute_jump!(result, rxn3a[1])
      @test result == expected

      result   = copy(x)
      expected = x + [-3, -1, -2, 1, 0]
      execute_jump!(result, rxn3b[1])
      @test result == expected
    end

    @testset "rate" begin
      @test rate(rxn0[1], x, params) ≈ params[1]
      @test rate(rxn0[2], x, params) ≈ params[1]

      @test rate(rxn1[1], x, params) ≈ x[1] * params[2]
      @test rate(rxn1[2], x, params) ≈ x[1] * params[2]
      @test rate(rxn1[3], x, params) ≈ x[1] * params[2]

      @test rate(rxn2a[1], x, params) ≈ x[1] * x[2] * params[3]
      @test rate(rxn2b[1], x, params) ≈ 1//2 * x[2] * (x[2] - 1) * params[4]

      @test rate(rxn3a[1], x, params) ≈ x[1] * x[2] * x[3] * params[1]
      @test rate(rxn3b[1], x, params) ≈ x[1] * (x[1] - 1) * (x[1] - 2) * x[2] * x[3] * (x[3] - 1) * params[1] / 12
    end
  end

  @testset "building from a Network" begin
    model = Network("silly test")

    # species definitions
    model <= Species("A", x[1])
    model <= Species("B", x[2])
    model <= Species("C", x[3])
    model <= Species("D", x[4])
    model <= Species("E", x[5])

    # order 0
    model <= Reaction("R1", params[1], "0 --> A")
    model <= Reaction("R2", params[1], "0 --> 4 * B")

    # order 1
    model <= Reaction("R3", params[2], "A --> 0")
    model <= Reaction("R4", params[2], "A --> B")
    model <= Reaction("R5", params[2], "A --> 2 * C + 2 * E")

    # order 2a
    model <= Reaction("R6", params[3], "A + B --> C")

    # order 2b
    model <= Reaction("R7", params[4], "2 * B --> A")

    # order 3
    model <= Reaction("R8", params[1], "A + B + C --> D")
    model <= Reaction("R9", params[1], "3 * A + B + 2 * C --> D")

    m = ReactionSystem(model)

    @testset "MassActionOrder0" begin
      expected = ReactionStruct(MassAction(), order0[1][1], order0[1][2], 1)

      @test m.reactions[1] == expected
      @test m.reactions[1].paramidx == expected.paramidx

      expected = ReactionStruct(MassAction(), order0[2][1], order0[2][2], 2)
      @test m.reactions[2] == expected
      @test m.reactions[2].paramidx == expected.paramidx
    end

    @testset "MassActionOrder1" begin
      expected = ReactionStruct(MassAction(), order1[1][1], order1[1][2], 3)

      @test m.reactions[3] == expected
      @test m.reactions[3].paramidx == expected.paramidx

      expected = ReactionStruct(MassAction(), order1[2][1], order1[2][2], 4)

      @test m.reactions[4] == expected
      @test m.reactions[4].paramidx == expected.paramidx

      expected = ReactionStruct(MassAction(), order1[3][1], order1[3][2], 5)

      @test m.reactions[5] == expected
      @test m.reactions[5].paramidx == expected.paramidx
    end

    @testset "MassActionOrder2A" begin
      expected = ReactionStruct(MassAction(), order2a[1][1], order2a[1][2], 6)

      @test m.reactions[6] == expected
      @test m.reactions[6].paramidx == expected.paramidx
    end

    @testset "MassActionOrder2B" begin
      expected = ReactionStruct(MassAction(), order2b[1][1], order2b[1][2], 7)

      @test m.reactions[7] == expected
      @test m.reactions[7].paramidx == expected.paramidx
    end

    @testset "MassActionOrderN" begin
      expected = ReactionStruct(MassAction(), order3a[1][1], order3a[1][2], 8)

      @test m.reactions[8] == expected
      @test m.reactions[8].paramidx == expected.paramidx

      expected = ReactionStruct(MassAction(), order3b[1][1], order3b[1][2], 9)
      @test m.reactions[9] == expected
      @test m.reactions[9].paramidx == expected.paramidx
    end

    @testset "fire_reaction!" begin
      result = copy(x); expected = x + [1, 0, 0, 0, 0]
      execute_jump!(result, m, 1)
      @test result == expected

      result = copy(x); expected = x + [0, 4, 0, 0, 0]
      execute_jump!(result, m, 2)
      @test result == expected

      result = copy(x); expected = x + [-1, 0, 0, 0, 0]
      execute_jump!(result, m, 3)
      @test result == expected

      result = copy(x); expected = x + [-1, 1, 0, 0, 0]
      execute_jump!(result, m, 4)
      @test result == expected

      result = copy(x); expected = x + [-1, 0, 2, 0, 2]
      execute_jump!(result, m, 5)
      @test result == expected

      result = copy(x); expected = x + [-1, -1, 1, 0, 0]
      execute_jump!(result, m, 6)
      @test result == expected

      result = copy(x); expected = x + [1, -2, 0, 0, 0]
      execute_jump!(result, m, 7)
      @test result == expected

      result = copy(x); expected = x + [-1, -1, -1, 1, 0]
      execute_jump!(result, m, 8)
      @test result == expected

      result = copy(x); expected = x + [-3, -1, -2, 1, 0]
      execute_jump!(result, m, 9)
      @test result == expected
    end

    @testset "rate" begin
      @test rate(m, x, 1) ≈ params[1]
      @test rate(m, x, 2) ≈ params[1]

      @test rate(m, x, 3) ≈ x[1] * params[2]
      @test rate(m, x, 4) ≈ x[1] * params[2]
      @test rate(m, x, 5) ≈ x[1] * params[2]

      @test rate(m, x, 6) ≈ x[1] * x[2] * params[3]
      @test rate(m, x, 7) ≈ 0.5 * x[2] * (x[2] - 1) * params[4]

      @test rate(m, x, 8) ≈ x[1] * x[2] * x[3] * params[1]
      @test rate(m, x, 9) ≈ x[1] * (x[1] - 1) * (x[1] - 2) * x[2] * x[3] * (x[3] - 1) * params[1] / 12
    end

    @testset "rate_derivative" begin
      # Order 0:
      @test rate_derivative(m, x, 1, 1) == 0

      # Order 1: A --> B
      @test rate_derivative(m, x, 1, 4) ≈ params[2]
      @test rate_derivative(m, x, 2, 4) == 0

      # Order 2A: A + B --> C
      @test rate_derivative(m, x, 1, 6) ≈ x[2] * params[3]
      @test rate_derivative(m, x, 2, 6) ≈ x[1] * params[3]
      @test rate_derivative(m, x, 3, 6) == 0

      # Order 2B: 2*B --> A
      @test rate_derivative(m, x, 2, 7) ≈ (x[2] - 0.5) * params[4]
      @test rate_derivative(m, x, 3, 7) == 0

      # Order 3: A + B + C --> D
      @test rate_derivative(m, x, 1, 8) ≈ x[2] * x[3] * params[1]
      @test rate_derivative(m, x, 2, 8) ≈ x[1] * x[3] * params[1]
      @test rate_derivative(m, x, 3, 8) ≈ x[1] * x[2] * params[1]
      @test rate_derivative(m, x, 4, 8) == 0

      # Order N: 3 * A + B + 2 * C --> D
      @test rate_derivative(m, x, 1, 9) ≈ (0.5*x[1]^2 - x[1] + 1/3) * x[2] * x[3] * (x[3] - 1) * params[1] / 2
      @test rate_derivative(m, x, 2, 9) ≈ x[1] * (x[1] - 1) * (x[1] - 2) * x[3] * (x[3] - 1) * params[1] / 12
      @test rate_derivative(m, x, 3, 9) ≈ x[1] * (x[1] - 1) * (x[1] - 2) * x[2] * (x[3] - 0.5) * params[1] / 6
      @test rate_derivative(m, x, 4, 9) == 0
    end
  end
end
