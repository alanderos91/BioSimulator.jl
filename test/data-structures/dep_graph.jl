@testset "DependencyGraph" begin
  model1 = kendall()
  model1_rxnrxn = [
    [1, 2], # X --> X + X
    [1, 2], # X --> 0
    [1, 2]  # 0 --> X
  ]
  model1_spcrxn = [
    [1, 2]
  ]

  model2 = neg_autoreg() # gene P2_gene RNA P P2
  model2_rxnrxn = [
    [1, 2, 3, 6], # 1. gene + P2 --> P2_gene
    [1, 2, 3, 6], # 2. P2_gene --> gene + P2
    [4, 7],       # 3. gene --> gene + RNA  / gene does not change number
    [5, 8],       # 4. RNA --> RNA + P      / RNA does not change number
    [1, 5, 6, 8], # 5. P + P --> P2
    [1, 5, 6, 8], # 6. P2 --> P + P
    [4, 7],       # 7. RNA --> 0
    [5, 8]        # 8. P --> 0
  ]
  model2_spcrxn = [
    [1, 3],
    [2],
    [4, 7],
    [5, 8],
    [1, 6]
  ]

  @testset "DGVector{$(T)}" for T in [DGView, DGLazy]
    @testset "Birth-Death-Immigration Process" begin
      @testset "reaction-reaction" begin
        dg = rxnrxn_depgraph(T(), model1)
        
        for j in eachindex(model1_rxnrxn)
          @test collect(dependents(dg, j)) == model1_rxnrxn[j]
        end
      end

      @testset "species-reaction" begin
        dg = spcrxn_depgraph(T(), model1)

        for i in eachindex(model1_spcrxn)
          @test collect(dependents(dg, i)) == model1_spcrxn[i]
        end
      end
    end

    @testset "Negative Gene Autoregulation" begin
      @testset "reaction-reaction" begin
        dg = rxnrxn_depgraph(T(), model2)

        for j in eachindex(model2_rxnrxn)
          @test collect(dependents(dg, j)) == model2_rxnrxn[j]
        end
      end

      @testset "species-reaction" begin
        dg = spcrxn_depgraph(T(), model2)

        for i in eachindex(model2_spcrxn)
          @test collect(dependents(dg, i)) == model2_spcrxn[i]
        end
      end
    end
  end
end