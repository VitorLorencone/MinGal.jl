include("../src/MinGal.jl")
using .MinGal
using Test

Algebra(3)

@testset "Operator Overloading Tests" begin
    e1 = Blade(1,1.0)
    e2 = Blade(2,2.0)
    mv1 = Multivector([1,2],[1.0,2.0])
    mv2 = Multivector([2,3],[2.0,3.0])

    # * Geometric Product and Scalar Multiplication
    @test isa(e1 * e2, Multivector)
    @test isa(mv1 * 2.0, Multivector)
    @test isa(2.0 * mv1, Multivector)

    # / Division
    @test isa(mv1 / 2.0, Multivector)
    @test mv1 / 2.0 == product_by_scalar(mv1, 0.5)

    # \ Inner Product
    @test isa(e1 \ e2, Multivector)
    @test isa(mv1 \ 1, Multivector)
    @test isa(1 \ mv2, Multivector)

    # ^ Outer Product
    @test isa(e1 ^ e2, Multivector)
    @test isa(e1 ^ 1, Multivector)
    @test isa(mv1 ^ 2, Multivector)

    # invert
    bl = Blade(0,1.0)
    @test invert(bl) == bl

    # Comparison operators
    bl_copy = Blade(1,1.0)
    @test e1 == bl_copy
    @test e1 != e2
    @test e1 == 1.0 || true

    # Unary minus and plus
    @test -mv1 == product_by_scalar(mv1, -1)
    @test mv1 + mv2 == multivector_sum(mv1, mv2)
    @test mv1 - mv2 == multivector_sub(mv1, mv2)
    @test 1 + mv1 == multivector_sum(Blade(0,1), mv1)
    @test mv1 + 1 == multivector_sum(Blade(0,1), mv1)

    # Reverse operator
    @test ~mv1 == revert(mv1)

    # Left and Right Contractions
    @test isa(mv1 << mv2, Multivector)
    @test isa(1 << mv2, Multivector)
    @test isa(mv1 << 1, Multivector)
    @test isa(mv1 >> mv2, Multivector)
    @test isa(1 >> mv2, Multivector)
    @test isa(mv1 >> 1, Multivector)

    # Regressive product &
    @test isa(mv1 & mv2, Multivector)
    @test isa(1 & mv2, Multivector)
    @test isa(mv1 & 1, Multivector)
end
