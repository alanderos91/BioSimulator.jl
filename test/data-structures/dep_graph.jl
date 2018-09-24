@testset "DependencyGraph" begin
  model1 = kendall()
  model1_deps = [
    [1, 2], # X --> X + X
    [1, 2], # X --> 0
    [1, 2]  # 0 --> X
  ]

  model2 = neg_autoreg()
  model2_deps = [
    [1, 2, 3, 6], # gene + P2 --> P2_gene
    [1, 2, 3, 6], # gene + P2 --> P2_gene
    [4, 7],       # gene --> gene + RNA  / gene does not change number
    [5, 8],       # RNA --> RNA + P      / RNA does not change number
    [1, 5, 6, 8], # P + P --> P2
    [1, 5, 6, 8], # P2 --> P + P
    [4, 7],       # RNA --> 0
    [5, 8]        # P --> 0
  ]

  @testset "DGVector{$(T)}" for T in [DGView, DGLazy]
    @testset "Birth-Death-Immigration Process" begin
      dg = DGVector{T}(model1)
      
      for j in eachindex(model1_deps)
        @test collect(dependents(dg, j)) == model1_deps[j]
      end
    end

    @testset "Negative Gene Autoregulation" begin
      dg = DGVector{T}(model2)

      for j in eachindex(model2_deps)
        @test collect(dependents(dg, j)) == model2_deps[j]
      end
    end
  end
end