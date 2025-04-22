include("../src/Mingal.jl")
using .Mingal
using Test

Complex = Algebra(0, 1, ["i"])
@test i*i == -1*id
@test 5*i*i*i == -5*i
@test (1+i)*(1+i) == 2*i
@test (1+i)*(1-i) == 2*id

R2 = Algebra(2,0,["i","j"])
@test i^i == 0*id
@test i\j == 0*id
@test i^j == ij
@test (i + 2*j)\(5*i+2*j) == 9*id

### Problems in changing algebra
### WARNING: redefinition of constant i. This may fail, cause incorrect answers, or produce other errors.