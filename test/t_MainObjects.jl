include("../src/MinGal.jl")
using .MinGal
using Test
using SparseArrays

Algebra(3)

@testset "Blade Tests" begin

    bl1 = Blade(0, 2.5)
    @test isa(bl1, Blade)
    @test bitmap(bl1) == 0
    @test scalar(bl1) == 2.5

    bl2 = Blade(7)
    @test isa(bl2, Blade)
    @test scalar(bl2) == 7

    mv = Multivector([0], [3.5])
    bl3 = Blade(mv)
    @test isa(bl3, Blade)
    @test scalar(bl3) == 3.5
    @test bitmap(bl3) == 0

    bl_copy = Blade(bl1)
    @test bl_copy === bl1
end

@testset "Multivector Tests" begin
    mv1 = Multivector([0, 1, 3], [1.0, -2.0, 3.5])
    @test isa(mv1, Multivector)
    @test length(mv1.blade_array.nzind) == 3

    bl = Blade(2, 4.0)
    mv2 = Multivector(bl)
    @test isa(mv2, Multivector)
    @test scalar(Blade(mv2)) == 4.0

    bl_list = [Blade(0, 1.0), Blade(1, -1.0), Blade(3, 2.5)]
    mv3 = Multivector(bl_list)
    @test isa(mv3, Multivector)
    @test length(mv3.blade_array.nzind) == 3
end

@testset "Exception Tests" begin
    mv_multi = Multivector([0,1], [1.0,2.0])
    @test_throws DomainError Blade(mv_multi)
end