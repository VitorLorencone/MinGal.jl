include("../src/MinGal.jl")
using .MinGal
using Test

# AlgebraFull
@testset "AlgebraFull Tests" begin
    Al3D = MinGal.create_algebra(3)
    
    @test isa(Al3D, MinGal.AlgebraFull)
    @test Al3D.p == 3
    @test Al3D.q == 0
    @test Al3D.r == 0
    @test length(Al3D.symbols) == 3
    @test length(Al3D.basis_bit_order) == 8
    @test Al3D.max == 8
    
    @test_throws DomainError MinGal.create_algebra(-1)
    @test_throws DomainError MinGal.create_algebra(3, -2)
    @test_throws DomainError MinGal.create_algebra(3, 0, -1)
    
    custom_symbols = ["x", "y", "z"]
    AlCustom = MinGal.create_algebra(3, 0, 0, custom_symbols)
    @test AlCustom.symbols == custom_symbols
    @test length(AlCustom.basis) == 8
end

# AlgebraMin
@testset "AlgebraMin Tests" begin
    AlMin = MinGal.create_algebra_min(2, 1)
    
    @test isa(AlMin, MinGal.AlgebraMin)
    @test AlMin.p == 2
    @test AlMin.q == 1
    @test AlMin.r == 0
    @test length(AlMin.metric) == 3
    @test AlMin.metric == Int8[1, 1, -1]
    @test AlMin.max == 1<<3
    
    @test_throws DomainError MinGal.create_algebra_min(-1)
    @test_throws DomainError MinGal.create_algebra_min(2, -3)
    @test_throws DomainError MinGal.create_algebra_min(2, 0, -5)
end