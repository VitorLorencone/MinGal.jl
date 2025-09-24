include("Algebra.jl")
import SparseArrays
using .SparseArrays

abstract type GAType end

"""
    Blade(blade_array)

Struct that creates the Blade object. The bitmap that specifies what basis vectors are present in this blade 
and The scalar of the basis blade.

# Arguments
- `blade_array::SparseArrays.SparseVector{Number, Int}` : An sparse vector with the internal values of 
basis blades and their scalars.


"""
struct Blade <: GAType
    blade_array::SparseArrays.SparseVector{Number, Int}
end

"""
    Multivector(blade_array)

Struct that creates the multivector object.

# Arguments
- `blade_array::SparseArrays.SparseVector{Number, Int}` : An sparse vector with the internal values of 
basis blades and their scalars.

"""
struct Multivector <: GAType
    blade_array::SparseArrays.SparseVector{Number, Int}
end

function Base.length(mv::GAType)::Int
    return length(mv.blade_array.nzind)
end

"""
    Blade(bitmap, scalar)::Blade ||
    Blade(mv::Multivector)::Blade

Creates a Blade based on bitmap and Scalars or Converts a multivector into a blade. It may throw an error

# Arguments
- `bitmap::Int` : The bitmap that specifies what basis vectors are present in this blade 
- `scalar::Number` : The scalar of the basis blade.
- `mv::Multivector`

# Return
Returns a Blade.

"""
function Blade(bitmap::Int, scalar)
    values = sparsevec([bitmap+1], [scalar], gb_current_algebra.max)
    return Blade(values)
end

function Blade(mv::Multivector)::Blade
    if length(mv) != 1
        throw(DomainError(mv, "You cannot convert this multivector to a blade."))
    end
    return Blade(mv.blade_array.nzind[1]-1, mv.blade_array.nzval[1])
end

"""
    bitmap(bl::Blade)::Int

Returns the bitmap of the Blade

# Arguments
- `bl::Blade`

# Return
The bitmap

"""

function bitmap(bl::Blade)::Int
    return bl.blade_array.nzind[1]-1
end

"""
    scalar(bl::Blade)::Number

Returns the scalar of the Blade

# Arguments
- `bl::Blade`

# Return
The scalar

"""

function scalar(bl::Blade)::Number
    return bl.blade_array.nzval[1]
end

"""
    Multivector(base_vectors, scalars) ||
    Multivector(bl::Blade) ||
    Multivector(blades::Array{Blade})

Constructor function for creating multivectors.

# Arguments
- `base_vectors::Array{Int}` : An array of integers, representing the actual basis blade that exists in this object in bit order.
- `scalars::Array{Number}` : An array of Floats, representing the scalars of each basis blade in bit order.
||
- `bl::Blade` : A Blade, to convert to Multivector.
||
- `blades::Array{Blade}` : An array of Blades, to convert to Multivector.

# Return
Returns a Multivector.

"""
function Multivector(base_vectors::Array, scalars::Array)::Multivector
    values = sparsevec(base_vectors .+ 1, scalars, gb_current_algebra.max)
    return Multivector(values)
end

function Multivector(bl::Blade)::Multivector
    values = sparsevec([bitmap(bl)+1], [scalar(bl)], gb_current_algebra.max)
    return Multivector(values)
end

function Multivector(blades::Array{Blade})::Multivector
    bitmaps = []
    scalars = []
    for bl in blades
        push!(bitmaps, bitmap(bl))
        push!(scalars, scalar(bl))
    end
    values = sparsevec(bitmaps, scalars, gb_current_algebra.max)
    return Multivector(values)
end

function str_bl(n::Int, al::AlgebraFull)::String
    return gb_current_algebra.basis_bit_order[n]
end

function str_bl(n::Int, al::AlgebraMin)::String
    n -= 1
    ans = ""
    for i in 1:sizeof(n) * 8
        if (n >> (i - 1)) & 1 == 1
            ans *= gb_current_algebra.symbols[i]
        end
    end
    return ans
end

# Pretty print a GAType
function Base.show(io::IO, mv::GAType)
    bitmaps = mv.blade_array.nzind
    scalars = mv.blade_array.nzval
    ans = ""

    for i in eachindex(bitmaps)
        if bitmaps[i] > 1 && scalars[i] > 0
            ans = ans * "+ " * string(round(scalars[i], digits=3)) * "*" * str_bl(bitmaps[i], gb_current_algebra) * " "
        elseif bitmaps[i] > 1 && scalars[i] < 0
            ans = ans * "- " * string(round(abs(scalars[i]), digits=3)) * "*" * str_bl(bitmaps[i], gb_current_algebra) * " "
        elseif bitmaps[i] == 1 && scalars[i] > 0
            ans = ans * string(scalars[i]) * " "
        elseif bitmaps[i] == 1 && scalars[i] < 0
            ans = ans * "- " * string(abs(scalars[i])) * " "
        elseif bitmaps[i] == 1 && scalars[i] == 0 && length(bitmaps) == 1
            ans = "0.0 "
        end
    end
    if length(ans) > 0
        if ans[1] == '+'
            ans = ans[3:end]
        end
    else
        ans = "0.0 "
    end
    ans = ans[1:end-1]
    print(ans)
end