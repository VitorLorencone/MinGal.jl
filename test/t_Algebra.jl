include("../src/Mingal.jl")
using .Mingal
using Test

Al3D = Algebra(3)
@test Al3D.p == 3
@test Al3D.q == 0
@test Al3D.VectorBasis == ["e1", "e2", "e3"]
@test Mingal.CurrentAlgebra.Basis == [("1", 1), ("e1", 2), ("e2", 3), ("e3", 4), ("e1e2", 5), ("e1e3", 6), ("e2e3", 7), ("e1e2e3", 8)]
@test Mingal.CurrentAlgebra.Indexes == [[0], [1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]]
@test_throws "The parameter p must be greater than 0" Algebra(-1)
@test_throws "The parameter q must be greater than 0" Algebra(3,-5)