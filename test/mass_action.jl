import BioSimulator: DenseReactionSystem,
                     SparseReactionSystem,
                     propensities,
                     scaled_rates,
                     update_all_propensities!,
                     PVec

∂ = BioSimulator.compute_mass_action_deriv

function mass_action_tests(Xt, r)
  X  = Xt[1]
  Y  = Xt[2]
  XY = Xt[3]

  k   = scaled_rates(r)
  k0  = k[1]
  k1  = k[2]
  k2a = k[3]
  k2b = k[4]

  a = propensities(r)

  Xt0 = zero(Xt)
  update_all_propensities!(r, Xt0)

  print("  Zero-Order:             ")
  @test a[1] == k0
  @test ∂(Xt0, r, 1, 1) == 0.0
  @test ∂(Xt0, r, 1, 2) == 0.0
  @test ∂(Xt0, r, 1, 3) == 0.0
  println("Passed")

  print("  First-Order:            ")
  @test a[2] == 0.0
  @test ∂(Xt0, r, 2, 1) == k1
  @test ∂(Xt0, r, 2, 2) == 0.0
  @test ∂(Xt0, r, 2, 3) == 0.0
  println("Passed")

  print("  Second-Order, Distinct: ")
  @test a[3] == 0.0
  @test ∂(Xt0, r, 3, 1) == 0.0
  @test ∂(Xt0, r, 3, 2) == 0.0
  @test ∂(Xt0, r, 3, 3) == 0.0
  println("Passed")

  print("  Second-Order, Repeated: ")
  @test a[4] == 0.0
  @test ∂(Xt0, r, 4, 1) == -k2b
  @test ∂(Xt0, r, 4, 2) == 0.0
  @test ∂(Xt0, r, 4, 3) == 0.0
  println("Passed")

  update_all_propensities!(r, Xt)

  print("  Zero-Order:             ")
  @test a[1] == k0
  @test ∂(Xt, r, 1, 1) == 0.0
  @test ∂(Xt, r, 1, 2) == 0.0
  @test ∂(Xt, r, 1, 3) == 0.0
  println("Passed")

  print("  First-Order:            ")
  @test a[2] == k1 * X
  @test ∂(Xt, r, 2, 1) == k1
  @test ∂(Xt, r, 2, 2) == 0.0
  @test ∂(Xt, r, 2, 3) == 0.0
  println("Passed")

  print("  Second-Order, Distinct: ")
  @test a[3] == k2a * X * Y
  @test ∂(Xt, r, 3, 1) == k2a * Y
  @test ∂(Xt, r, 3, 2) == k2a * X
  @test ∂(Xt, r, 3, 3) == 0.0
  println("Passed")

  print("  Second-Order, Repeated: ")
  @test a[4] == k2b * X * (X - 1)
  @test ∂(Xt, r, 4, 1) == k2b * (2 * X - 1)
  @test ∂(Xt, r, 4, 2) == 0.0
  @test ∂(Xt, r, 4, 3) == 0.0
  println("Passed")
end

##### Setup #####
Xt = [100, 500, 1] # (X, Y, XY)

pre = transpose([
  0 0 0; # Zero Order
  1 0 0; # First Order
  1 1 0; # Second Order, Type A
  2 0 0  # Second Order, Type B
])

post = transpose([
  1 0 0; # Zero Order
  0 0 0; # First Order
  0 0 1; # Second Order, Type A
  0 1 0  # Second Order, Type B
])

V = post - pre
U = pre

k = [0.1, 0.01, 0.25, 1/2 * 0.5]
a = PVec{Float64}(4)
dg = Array{Int}[] # unused

##### DenseReactionSystem #####
println("- DenseReactionSystem")
mass_action_tests(Xt, DenseReactionSystem(V, U, k, a, dg))

##### SparseReactionSystem #####
println("- SparseReactionSystem")
mass_action_tests(Xt, SparseReactionSystem(V, U, k, a, dg))
