import BioSimulator: DenseReactionSystem,
                     SparseReactionSystem,
                     propensities,
                     compute_propensities!

∂ = BioSimulator.mass_action_deriv

function mass_action_tests(Xt, r, parameters)
  X  = Xt[1]
  Y  = Xt[2]
  XY = Xt[3]

  k0  = value(parameters[:k0])
  k1  = value(parameters[:k1])
  k2a = value(parameters[:k2a])
  k2b = value(parameters[:k2b])

  a = propensities(r)

  Xt0 = zero(Xt)
  compute_propensities!(r, Xt0, parameters)

  print("  Zero-Order:             ")
  @test a[1] == k0
  @test ∂(Xt0, r, parameters, 1, 1) == 0.0
  @test ∂(Xt0, r, parameters, 1, 2) == 0.0
  @test ∂(Xt0, r, parameters, 1, 3) == 0.0
  println("Passed")

  print("  First-Order:            ")
  @test a[2] == 0.0
  @test ∂(Xt0, r, parameters, 2, 1) == k1
  @test ∂(Xt0, r, parameters, 2, 2) == 0.0
  @test ∂(Xt0, r, parameters, 2, 3) == 0.0
  println("Passed")

  print("  Second-Order, Distinct: ")
  @test a[3] == 0.0
  @test ∂(Xt0, r, parameters, 3, 1) == 0.0
  @test ∂(Xt0, r, parameters, 3, 2) == 0.0
  @test ∂(Xt0, r, parameters, 3, 3) == 0.0
  println("Passed")

  print("  Second-Order, Repeated: ")
  @test a[4] == 0.0
  @test ∂(Xt0, r, parameters, 4, 1) == -k2b
  @test ∂(Xt0, r, parameters, 4, 2) == 0.0
  @test ∂(Xt0, r, parameters, 4, 3) == 0.0
  println("Passed")

  compute_propensities!(r, Xt, parameters)

  print("  Zero-Order:             ")
  @test a[1] == k0
  @test ∂(Xt, r, parameters, 1, 1) == 0.0
  @test ∂(Xt, r, parameters, 1, 2) == 0.0
  @test ∂(Xt, r, parameters, 1, 3) == 0.0
  println("Passed")

  print("  First-Order:            ")
  @test a[2] == k1 * X
  @test ∂(Xt, r, parameters, 2, 1) == k1
  @test ∂(Xt, r, parameters, 2, 2) == 0.0
  @test ∂(Xt, r, parameters, 2, 3) == 0.0
  println("Passed")

  print("  Second-Order, Distinct: ")
  @test a[3] == k2a * X * Y
  @test ∂(Xt, r, parameters, 3, 1) == k2a * Y
  @test ∂(Xt, r, parameters, 3, 2) == k2a * X
  @test ∂(Xt, r, parameters, 3, 3) == 0.0
  println("Passed")

  print("  Second-Order, Repeated: ")
  @test a[4] == k2b * X * (X - 1)
  @test ∂(Xt, r, parameters, 4, 1) == k2b * (2 * X - 1)
  @test ∂(Xt, r, parameters, 4, 2) == 0.0
  @test ∂(Xt, r, parameters, 4, 3) == 0.0
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

v = post - pre
u = pre

k_id = [:k0, :k1, :k2a, :k2b]
k    = [0.1, 0.01, 0.25, 0.5]

parameters = Dict{Symbol,Parameter}(
  k_id[1] => Parameter(k_id[1], k[1]),
  k_id[2] => Parameter(k_id[2], k[2]),
  k_id[3] => Parameter(k_id[3], k[3]),
  k_id[4] => Parameter(k_id[4], k[4])
)

##### DenseReactionSystem #####
println("- DenseReactionSystem")
mass_action_tests(Xt, DenseReactionSystem(v, u, k_id), parameters)

##### SparseReactionSystem #####
println("- SparseReactionSystem")
mass_action_tests(Xt, SparseReactionSystem(v, u, k_id), parameters)
