include("../src/MinGal.jl")
using .MinGal
using Test

Al3D = Algebra(3)

@testset "Canonical Reordering Sign Tests" begin
    bl1 = Blade(1,1.0)
    bl2 = Blade(2,1.0)
    bl3 = Blade(3,1.0)

    s1 = MinGal.canonical_reordering_sign(bl1, bl2)
    s2 = MinGal.canonical_reordering_sign(bl2, bl1)

    @test s1 == 1.0 || s1 == -1.0
    @test s2 == 1.0 || s2 == -1.0
end

@testset "Geometric, Outer, Inner Products" begin
    e1 = Blade(1,1.0)
    e2 = Blade(2,1.0)

    gp = geometric_product(e1,e2)
    @test isa(gp, Multivector)

    op = outer_product(e1,e2)
    @test isa(op, Multivector)
    @test scalar(Blade(op)) == 1.0

    ip = inner_product(e1,e2)
    @test isa(ip, Multivector)
end

@testset "Sum, Sub, Scalar Multiplication" begin
    e1 = Blade(1,2.0)
    e2 = Blade(2,3.0)
    mv1 = Multivector([1], [2.0])
    mv2 = Multivector([2], [3.0])

    sum_mv = multivector_sum(mv1, mv2)
    @test length(sum_mv) == 2

    sub_mv = multivector_sub(mv1, mv2)
    @test length(sub_mv) == 2

    prod_mv = product_by_scalar(mv1, 5.0)
    @test scalar(Blade(prod_mv)) == 10.0
end

@testset "Reverse and Grade Involution" begin
    e1 = Blade(1,1.0)
    e2 = Blade(3,1.0)
    mv = Multivector([1,3],[1.0,2.0])

    rev_mv = revert(mv)
    inv_mv = grade_involution(mv)

    @test isa(rev_mv, Multivector)
    @test isa(inv_mv, Multivector)
end

@testset "Left and Right Contractions" begin
    e1 = Blade(1,1.0)
    e2 = Blade(3,1.0)
    lc = left_contraction(e1,e2)
    rc = right_contraction(e2,e1)

    @test isa(lc, Multivector)
    @test isa(rc, Multivector)
end

@testset "Grade Selection, Dual, Undual" begin
    e1 = Blade(1,1.0)
    mv = Multivector([1,3],[1.0,2.0])

    gs = grade_selection(mv,1)
    @test length(gs) == 1

    dual_mv = dual(mv)
    undual_mv = undual(dual_mv)
    @test isa(dual_mv, Multivector)
    @test isa(undual_mv, Multivector)
end

@testset "Regressive Product" begin
    e1 = Blade(1,1.0)
    e2 = Blade(2,1.0)
    rp = regressive_product(e1,e2)
    @test isa(rp, Multivector)
end

@testset "Scalar / Blade Checks" begin
    e1 = Blade(1,1.0)
    mv = Multivector([0],[2.0])

    @test isScalar(mv) == true
    @test isBlade(e1) == true
end

@testset "Exponential of GAType" begin
    e1 = Blade(1,1.0)
    exp_mv = exp_ga(e1)
    @test isa(exp_mv, Multivector)
end

@testset "Mix" begin
    @test grade(e3) == 1
    @test grade(e1e3) == 2
    @test grade(id) == 0
    @test e3*e1 == -e1e3
    @test -2*e1*5.3 == -10.6*e1
    @test e1e2\e1 == -e2
    @test e1^e1e2 == 0*id
    @test 5*e3*7*e1e3 == -35*e1
    @test e1+e2 == Multivector([1, 2],[1, 1])
    @test e1-e2 == Multivector([1, 2],[1, -1])
    @test e1+e2-e1 == e2
    @test e1+1-e1 == 1.0*id
    @test 2*e1+3*e1+4*e3-7*e3-6*e1e2e3+2.5*e1e2 == 5*e1 - 3*e3 - 6*e1e2e3 + 2.5*e1e2
    @test e1*(3+1.5*e3) == 3*e1 + 1.5*e1e3
    @test e1\(3+1.5*e2) == 3*e1
    @test e1^(3+1.5*e1) == 3*e1
    @test e1^(0+e1+e2) == e1e2
    @test (0+e1)^e1 == 0*id
    @test (1+e1+e2+e3+e1e2+e1e3+e2e3+e1e2e3)*(1+e1+e2+e3+e1e2+e1e3+e2e3+e1e2e3) == 4*e2 + 4*e1e2 + 4*e2e3 + 4*e1e2e3
    @test (1+e1+e2+e3+e1e2+e1e3+e2e3+e1e2e3)\(1+e1+e2+e3+e1e2+e1e3+e2e3+e1e2e3) == 4*e2 + 4*e1e2 + 4*e2e3 + 2*e1e2e3
    @test (1+e1+e2+e3+e1e2+e1e3+e2e3+e1e2e3)^(1+e1+e2+e3+e1e2+e1e3+e2e3+e1e2e3) == 1 + 2*e1 + 2*e2 + 2*e3 + 2*e1e2 + 2*e1e3 + 2*e2e3 + 4*e1e2e3
    @test e1e2 + 2*e3\((3*e1 + 4*e1e2) \ (e1e2e3)) == e1e2 - 6*e2 - 8
end

Algebra2 = Algebra(2,1,1)

@testset "GeometricFunctions Tests Algebra(2,1,1)" begin

    mv1 = Multivector([1,2],[1.0,2.0])
    mv2 = Multivector([2,3],[3.0,4.0])

    @test isa(e1 * e2, Multivector)
    @test isa(e1 ^ e2, Multivector)
    @test isa(e1 \ e3, Multivector)

    @test isa(mv1 * 2.0, Multivector)
    @test isa(2.0 * mv1, Multivector)
    @test isa(mv1 / 2.0, Multivector)

    e1copy = Blade(1,1.0)
    @test e0 == e1copy
    @test e1 != e2

    @test -mv1 == product_by_scalar(mv1,-1)
    @test mv1 + mv2 == multivector_sum(mv1,mv2)
    @test mv1 - mv2 == multivector_sub(mv1,mv2)

    @test ~mv1 == revert(mv1)

    @test isa(mv1 << mv2, Multivector)
    @test isa(mv1 >> mv2, Multivector)
    @test isa(e1 << e3, Multivector)
    @test isa(e2 >> e0, Multivector)

    @test isa(mv1 & mv2, Multivector)
    @test isa(e1 & e2, Multivector)

    dual_mv1 = dual(mv1)
    undual_mv1 = undual(dual_mv1)
    @test isa(dual_mv1, Multivector)
    @test isa(undual_mv1, Multivector)

    @test grade_projection(e1, 1) == e1
    @test grade_projection(e1, 2) == Multivector([0],[0])
end