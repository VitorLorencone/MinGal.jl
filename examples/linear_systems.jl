# None of them gets past D = 12, maybe with multiplicative
# implementation instead of the additive one???

using MinGal

"""
    solve_linear_systems(coeff::Array, rhs::Array)::Array

Function that solves a linear systems of equations using Geometric Algebra
with the Grassmann method

# Arguments
- `coeff::Array` : The matrix of coeficients, the left side of the equations.
- `rhs::Array` : The matrix of values, the right side of the equations.

# Return
An Array with the solutions.

"""
function solve_linear_systems(coeff::Array, rhs::Array)
    # assume the inputs are always correct
    D = size(coeff,1)
    Algebra(D)

    # Calculates the grassmann coeficients such as A1 = a1*e1 + a2*e2 +...+ an*en
    coeff_grass::GAVector = [sum(Blade(1<<(j-1), coeff[j, i]) for j in 1:D) for i in 1:D]

    # The same thing but for the right side T = b1*e1 + b2*e2 +...+ bn*en
    T::GAType = sum(Blade(1<<(j-1), rhs[j]) for j in 1:D)

    # The outer product of every element in coeff_grass, I = A1^A2^...^An
    I = invert(reduce(^, coeff_grass; init=Blade(1)))

    # for each variable x in the linear system, compute I^-1 * (A1^A2^...Ax-1^T^Ax+1^...^An)
    ans = zeros(D)
    for i in 1:D
        substituted_vecs = [j == i ? T : coeff_grass[j] for j in 1:D]
        partial = reduce(^, substituted_vecs; init=Blade(1))
        ans[i] = (I * partial)[0]
    end

    return ans
end

function solve_linear_systems_homogeneus(coeff::Array, rhs::Array)
    # assume the inputs are always correct
    D = size(coeff,1)
    Algebra(D, 0, 1)

    # Calculates the grassmann coeficients such as A1 = a*e1 + b*e2 +...- B1*e0
    coeff_grass = Multivector[]
    for i in 1:D
        mv = sum(Blade(1<<(j), coeff[i,j]) for j in 1:D) - rhs[i]*e0
        push!(coeff_grass, mv)
    end

    mv = dual(reduce(^, coeff_grass; init=Blade(1)))
    # mv = reduce(&, dual.(coeff_grass); init=eI) -> Alternative version, but worst

    if(mv == 0)
        error("Compatible and Indeterminate System - Infinite Solutions")
    elseif(mv[1] == 0)
        error("Impossible System - No Solution")
    end

    mv = mv/mv[1]
    ans = [mv[1<<k] for k in 1:D]
end

# ---------- Example 1 (3x3) ----------

A1 = [1 2 3;
     0 1 4;
     5 6 0]

B1 = [14, 14, 17]

# Output = [1.0, 2.0, 3.0]
println(solve_linear_systems_homogeneus(A1, B1))

# ---------- Example 2 (7x7) ----------

A2 = [
    1 0 2 1 0 3 1;
    0 1 4 0 1 2 0;
    2 1 0 3 1 0 1;
    1 2 1 1 0 1 2;
    0 1 2 0 1 0 3;
    3 0 1 2 0 1 1;
    1 1 0 1 2 3 0
]

B2 = [32, 25, 20, 19, 21, 20, 21]

# Output = [-2.394202898550725, -1.7768115942028986, 3.7594202898550724, 6.83768115942029, 
#           1.4492753623188406,  5.144927536231884,  4.602898550724638]
println(solve_linear_systems_homogeneus(A2, B2))

# ---------- Example 3 (10x10) ----------

D = 10
x_sol = collect(1:D)
A3 = rand(1:D, D, D)
B3 = A3 * x_sol

# Output = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
println(solve_linear_systems_homogeneus(A3, B3))