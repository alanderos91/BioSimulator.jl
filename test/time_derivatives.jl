import BioSimulator: DenseReactionSystem,
                     SparseReactionSystem,
                     propensities,
                     mean_derivatives!,
                     time_derivatives!

∂ = BioSimulator.compute_mass_action_deriv

function time_derivatives_tests(Xt, r)
  a = propensities(r)
  update_all_propensities!(a, r, Xt)

  dxdt = zeros(Float64, length(Xt))
  drdt = zeros(Float64, length(a))

  mean_derivatives!(dxdt, r)
  time_derivatives!(drdt, Xt, r, dxdt)

  print("  Species Mean Derivatives: ")
  val = a[1] * 1 + a[2] * -1 + a[3] * -1 + a[4] * -2
  @test dxdt[1] == val

  val = a[1] * 0 + a[2] * 0 + a[3] * -1 + a[4] * 1
  @test dxdt[2] == val

  val = a[1] * 0 + a[2] * 0 + a[3] * 1 + a[4] * 0
  @test dxdt[3] == val

  println("Passed")

  print("  Reaction Derivatives:     ")
  val = ∂(Xt, r, 1, 1) * dxdt[1] + ∂(Xt, r, 1, 2) * dxdt[2] + ∂(Xt, r, 1, 3) * dxdt[3]
  @test drdt[1] == val

  val = ∂(Xt, r, 2, 1) * dxdt[1] + ∂(Xt, r, 2, 2) * dxdt[2] + ∂(Xt, r, 2, 3) * dxdt[3]
  @test drdt[2] == val

  val = ∂(Xt, r, 3, 1) * dxdt[1] + ∂(Xt, r, 3, 2) * dxdt[2] + ∂(Xt, r, 3, 3) * dxdt[3]
  @test drdt[3] == val

  val = ∂(Xt, r, 4, 1) * dxdt[1] + ∂(Xt, r, 4, 2) * dxdt[2] + ∂(Xt, r, 4, 3) * dxdt[3]
  @test drdt[4] == val

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

k    = [0.1, 0.01, 0.25, 0.5]

##### DenseReactionSystem #####
println("- DenseReactionSystem")
time_derivatives_tests(Xt, DenseReactionSystem(V, U, k, a, dg))

##### SparseReactionSystem #####
println("- SparseReactionSystem")
time_derivatives_tests(Xt, SparseReactionSystem(V, U, k, a, dg))
