import BioSimulator: DenseReactionSystem,
                     SparseReactionSystem,
                     propensities,
                     compute_propensities!,
                     mean_derivatives!,
                     time_derivatives!

∂ = BioSimulator.mass_action_deriv

function time_derivatives_tests(Xt, r, parameters)
  a = propensities(r)
  compute_propensities!(r, Xt, parameters)

  dxdt = zeros(Float64, length(Xt))
  drdt = zeros(Float64, length(a))

  mean_derivatives!(dxdt, r)
  time_derivatives!(drdt, Xt, r, parameters, dxdt)

  print("  Species Mean Derivatives: ")
  val = a[1] * 1 + a[2] * -1 + a[3] * -1 + a[4] * -2
  @test dxdt[1] == val

  val = a[1] * 0 + a[2] * 0 + a[3] * -1 + a[4] * 1
  @test dxdt[2] == val

  val = a[1] * 0 + a[2] * 0 + a[3] * 1 + a[4] * 0
  @test dxdt[3] == val

  println("Passed")

  print("  Reaction Derivatives:     ")
  val = ∂(Xt, r, parameters, 1, 1) * dxdt[1] + ∂(Xt, r, parameters, 1, 2) * dxdt[2] + ∂(Xt, r, parameters, 1, 3) * dxdt[3]
  @test drdt[1] == val

  val = ∂(Xt, r, parameters, 2, 1) * dxdt[1] + ∂(Xt, r, parameters, 2, 2) * dxdt[2] + ∂(Xt, r, parameters, 2, 3) * dxdt[3]
  @test drdt[2] == val

  val = ∂(Xt, r, parameters, 3, 1) * dxdt[1] + ∂(Xt, r, parameters, 3, 2) * dxdt[2] + ∂(Xt, r, parameters, 3, 3) * dxdt[3]
  @test drdt[3] == val

  val = ∂(Xt, r, parameters, 4, 1) * dxdt[1] + ∂(Xt, r, parameters, 4, 2) * dxdt[2] + ∂(Xt, r, parameters, 4, 3) * dxdt[3]
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
time_derivatives_tests(Xt, DenseReactionSystem(v, u, k_id), parameters)

##### SparseReactionSystem #####
println("- SparseReactionSystem")
time_derivatives_tests(Xt, SparseReactionSystem(v, u, k_id), parameters)
