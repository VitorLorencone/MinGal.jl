include("../src/MinGal.jl")
using .MinGal
using Test

@testset "Complex Tests" begin
    Complex = Algebra(0,1,0,["i"])
    @test i*i == -1*id
    @test 5*i*i*i == -5*i
    @test (1+i)*(1+i) == 2*i
    @test (1+i)*(1-i) == 2*id
    @test i\i == -1
    @test i\id == i
end

@testset "R2 Tests" begin
    R2 = Algebra(2,0,0,["i","j"])
    @test i^i == 0*id
    @test i\j == 0*id
    @test i^j == ij
    @test (i + 2*j)\(5*i+2*j) == 9*id
end

@testset "Workspace Tests" begin

    AlFull = Algebra(2,1,1,nothing,"full")
    @test isa(AlFull, MinGal.AlgebraFull)
    @test length(AlFull.basis_bit_order) == 16

    AlMin = Algebra(2,1,1,nothing,"min")
    @test isa(AlMin, MinGal.AlgebraMin)
    @test length(AlMin.symbols) == 4

    AlSpecial = Algebra(2,1,1,nothing,"special")
    @test isa(AlSpecial, MinGal.AlgebraMin)

    @test isdefined(Main, :id)
    @test isdefined(Main, :eI)

    symbols_custom = ["", "a", "b", "c"]
    MinGal.create_symbols(symbols_custom)
    @test isdefined(Main, :a)
    @test isdefined(Main, :b)
    @test isdefined(Main, :c)
    @test isa(a, Multivector)
    @test isa(b, Multivector)
    @test isa(c, Multivector)

    MinGal.create_symbols_min(["x", "y", "z"])
    @test isdefined(Main, :x)
    @test isdefined(Main, :y)
    @test isdefined(Main, :z)
    @test isa(x, Multivector)
    @test isa(y, Multivector)
    @test isa(z, Multivector)

end

### Problems in changing algebra
### WARNING: redefinition of constant i. This may fail, cause incorrect answers, or produce other errors.