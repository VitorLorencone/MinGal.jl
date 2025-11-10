include("../src/MinGal.jl")
using .MinGal
using Test

Algebra(3)

@testset "Grade Tests" begin
    bl1 = Blade(0, 2.0)
    bl2 = Blade(3, 1.0) # bitmap 3 -> grade 2
    mv = Multivector([0,1,3], [1.0,2.0,3.0])

    @test grade(bl1) == 0
    @test grade(bl2) == 2
    @test grade(mv) == 2
end

@testset "Grade Projection Tests" begin
    bl = Blade(3, 1.5)
    proj1 = grade_projection(bl, 2)
    proj2 = grade_projection(bl, 1)

    @test isa(proj1, Multivector)
    @test scalar(Blade(proj1)) == 1.5
    @test scalar(Blade(proj2)) == 0.0

    mv = Multivector([3], [2.0])
    proj_mv = grade_projection(mv, 2)
    @test scalar(Blade(proj_mv)) == 2.0
end

@testset "Scalar Product Tests" begin
    bl1 = Blade(1, 1.0)
    bl2 = Blade(2, 1.0)
    bl3 = Blade(1, 2.0)

    @test MinGal.scalar_product(bl1, bl2) == 0

    @test MinGal.scalar_product(bl1, bl3) == 1

    bl_wrong = Blade(3, 1.0)
    @test_throws DomainError MinGal.scalar_product(bl_wrong, bl1)
end

@testset "Scalar Access Tests" begin
    mv = Multivector([0,1,3], [1.0,2.0,3.0])
    bl = Blade(1,1.0)

    # get_scalar
    @test MinGal.get_scalar(mv, 0) == 1.0
    @test MinGal.get_scalar(mv, bl) == 2.0

    # set_scalar
    MinGal.set_scalar(mv, 10.0, 0)
    @test MinGal.get_scalar(mv, 0) == 10.0
    MinGal.set_scalar(mv, 20.0, bl)
    @test MinGal.get_scalar(mv, bl) == 20.0

    # has_key
    @test MinGal.has_key(mv, 0) == true
    @test MinGal.has_key(mv, 2) == false
    @test MinGal.has_key(mv, bl) == true
end

@testset "Canonical Basis and Chain Tests" begin
    cb = canonical_basis()
    @test length(cb) == MinGal.gb_current_algebra.p + MinGal.gb_current_algebra.q + MinGal.gb_current_algebra.r
    coeffs = [1.0, 2.0, 3.0]
    chained = chain(coeffs)
    @test isa(chained, Blade) || isa(chained, Multivector)
end

@testset "Base function substitutions" begin
    mv = Multivector([0,1], [1.0,2.0])
    # getindex / setindex!
    @test mv[0] == 1.0
    mv[0] = 10.0
    @test mv[0] == 10.0

    # iterate
    iter = collect(mv)
    @test length(iter) == 2
    @test isa(iter[1], Blade)

    # keys / values
    @test keys(mv) == mv.blade_array.nzind .- 1
    @test values(mv) == mv.blade_array.nzval

    # haskey
    @test haskey(mv, 0) == true
end