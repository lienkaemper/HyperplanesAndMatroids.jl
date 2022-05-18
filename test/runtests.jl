using HyperplanesAndMatroids
using Test


@testset "HyperplanesAndMatroids.jl" begin
    # Write your tests here.
    h_degen = [1 0; -1 0; 0 1; -1 1; 1 1]
    M_non_uniform = HyperplanesAndMatroids.fromHyperplanes(h_degen)
    @testset "topes" begin
      @test hasproperty(M_non_uniform, :topes) == true
      @test issetequal(M_non_uniform.topes, HyperplanesAndMatroids.SignVector.(["+----", "+---+", "-+++-", "-++++", "+-+++", "+-+-+", "-+-+-", "-+---"]))
    end
    @testset "chirotope" begin
      @test M_non_uniform.chirotope == Dict(
        [1, 5] => 1,
        [2, 4] => -1,
        [4, 5] => -1,
        [1, 2] => 0,
        [1, 3] => 1,
        [2, 5] => -1,
        [3, 4] => 1,
        [1, 4] => 1,
        [2, 3] => -1,
        [3, 5] => -1)
      end
      @testset "cocircuits" begin
        @test issetequal(M_non_uniform.cocircuits[[3, 4, 5]], HyperplanesAndMatroids.SignVector.(["00---", "00+++"]))
        @test issetequal(M_non_uniform.cocircuits[[1,2,3,4]], HyperplanesAndMatroids.SignVector.(["+---0", "-+++0"]))
      end
      @testset "circuits" begin
        @test issetequal(M_non_uniform.circuits[[1,2]], HyperplanesAndMatroids.SignVector.(["++000", "--000"]))
        @test issetequal(M_non_uniform.circuits[[3,4,5]], HyperplanesAndMatroids.SignVector.(["00+--", "00-++"]))
      end
end
