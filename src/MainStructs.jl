include("Algebra.jl")
import SparseArrays
using .SparseArrays

# Abstract Geometric Type for Multivectors and Blades
abstract type AbstractGeometricAlgebraType end

"""
    Multivector(val)

Struct that creates the multivector object.

# Arguments
- `val::SparseArrays.SparseVector{Float64, Int64}` : An sparse vector with the internal values of basis blades and their scalars.

"""
mutable struct Multivector <: AbstractGeometricAlgebraType

    val::SparseArrays.SparseVector{Float64, Int64}

end

"""
    Blade(val)

Struct that creates the Blade object.

# Arguments
- `val::SparseArrays.SparseVector{Float64, Int64}` : An sparse vector with the internal values of basis blades and their scalars.

"""
mutable struct Blade <: AbstractGeometricAlgebraType

    val::SparseArrays.SparseVector{Float64, Int64}

end

"""
    Multivectors(baseVectors, scalars, Al)::AbstractGeometricAlgebraType

Constructor function for creating either blades or multivectors, done automatically.

# Arguments
- `baseVectors::Array` : An array of integers, representing the actual basis blade that exists in this object in order.
- `scalars::Array` : An array of integers, representing the scalars of each basis blade in order.
- `Al::Algebra` : The Algebra, it is setted as CurrentAlgebra.

# Return
Returns an AbstractGeometricAlgebraType.

"""
function Multivectors(baseVectors::Array, scalars::Array, Al::Algebra = CurrentAlgebra)

    max = 2^(Al.p + Al.q)

    for i in baseVectors
        # Internal Error
        @assert 0 <= i <= max
    end

    values = sparsevec(baseVectors, scalars, max)
    if length(baseVectors) == 1
        return Blade(values)
    else
        return Multivector(values)
    end
end

# Write an AbstractGeometricAlgebraType in REPL
function Base.show(io::IO, a::AbstractGeometricAlgebraType)
    ind = a.val.nzind
    val = a.val.nzval
    ans = ""

    for i in eachindex(ind)
        if ind[i] > 1 && val[i] > 0
            ans = ans * "+ " * string(round(val[i], digits=3)) * "*" * CurrentAlgebra.Basis[ind[i]][1] * " "
        elseif ind[i] > 1 && val[i] < 0
            ans = ans * "- " * string(round(abs(val[i]), digits=3)) * "*" * CurrentAlgebra.Basis[ind[i]][1] * " "
        elseif ind[i] == 1 && val[i] > 0
            ans = ans * string(val[i]) * " "
        elseif ind[i] == 1 && val[i] < 0
            ans = ans * "- " * string(val[i]) * " "
        elseif ind[i] == 1 && val[i] == 0 && length(ind) == 1
            ans = "0.0 "

        if ans[1] == '+'
            ans = ans[3:end]
        end

        ans = ans[1:end-1]

        print(ans)
end