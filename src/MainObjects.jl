include("Algebra.jl")
import SparseArrays
using .SparseArrays

abstract type GAType end
const GAVector = AbstractVector{<:GAType}
const GAArray{T<:GAType,N} = AbstractArray{T,N}

"""
    Blade(blade_array)

Struct that creates the Blade object. The bitmap that specifies what basis vectors are present in this blade 
and The scalar of the basis blade.

# Arguments
basis blades and their scalars.

"""
struct Blade <: GAType
    blade_array::SparseArrays.SparseVector{Number, Integer}
    algebra::Algebra
end

"""
    Multivector(blade_array)

Struct that creates the multivector object.

# Arguments
- `blade_array::SparseArrays.SparseVector{Number, Integer}` : An sparse vector with the internal values of 
basis blades and their scalars.

"""
struct Multivector <: GAType
    blade_array::SparseArrays.SparseVector{Number, Integer}
    algebra::Algebra
end

function Base.length(mv::GAType)::Integer
    dropzeros!(mv.blade_array)
    return length(mv.blade_array.nzind)
end

"""
    Blade(bitmap, scalar, al::Algebra)::Blade ||
    Blade(mv::Multivector)::Blade ||
    Blade(k::Number, al::Algebra)::Blade

Creates a Blade based on bitmap and Scalars or Converts a multivector into a blade. It may throw an error

# Arguments
- `bitmap::Integer` : The bitmap that specifies what basis vectors are present in this blade 
- `scalar::Number` : The scalar of the basis blade.
- `al::Algebra` : The current algebra on use.
||
- `mv::Multivector`
||
- `k::Number`
- `al::Algebra` : The current algebra on use.

# Return
Returns a Blade.

"""
function Blade(bitmap::Integer, scalar::Number, al::Algebra = gb_current_algebra)
    if typeof(bitmap) != typeof(al.max)
        bitmap = typeof(al.max)(bitmap)
    end
    
    values = sparsevec([bitmap+1], [scalar], al.max)
    return Blade(values, al)
end

function Blade(mv::Multivector)::Blade
    if length(mv) > 1
        throw(DomainError(mv, "You cannot convert this multivector to a blade."))
    elseif length(mv) == 0
        return Blade(0, mv.algebra)
    end
    return Blade(mv.blade_array.nzind[1]-1, mv.blade_array.nzval[1], mv.algebra)
end

function Blade(k::Number, al::Algebra = gb_current_algebra)::Blade
    bitmap = 1
    
    if typeof(al.max) == BigInt
        bitmap = typeof(al.max)(bitmap)
    end
    
    values = sparsevec([bitmap], [k], al.max)
    return Blade(values, al)
end

function Blade(bl::Blade)::Blade
    return bl
end

"""
    bitmap(bl::Blade)::Integer

Returns the bitmap of the Blade

# Arguments
- `bl::Blade`

# Return
The bitmap

"""

function bitmap(bl::Blade)::Integer
    return bl.blade_array.nzind[1]-1
end

function bitmap(mv::Multivector)::Integer
    return bitmap(Blade(mv))
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
    Multivector(base_vectors, scalars, al::Algebra) ||
    Multivector(bl::Blade) ||
    Multivector(blades::Array{Blade}, al::Algebra)

Constructor function for creating multivectors.

# Arguments
- `base_vectors::Vector{Integer}` : An array of integers, representing the actual basis blade that exists in this object in bit order.
- `scalars::Vector{Number}` : An array of Floats, representing the scalars of each basis blade in bit order.
- `al::Algebra` : The current algebra on use.
||
- `bl::Blade` : A Blade, to convert to Multivector.
||
- `blades::GAVector` : An array of Blades, to convert to Multivector.
- `al::Algebra` : The current algebra on use.

# Return
Returns a Multivector.

"""
function Multivector(base_vectors::Array, scalars::Array, al::Algebra = gb_current_algebra)::Multivector
    if typeof(al.max) == BigInt
        base_vectors = BigInt.(base_vectors)
    end
    
    values = sparsevec(base_vectors .+ 1, scalars, al.max)
    return Multivector(values, al)
end

function Multivector(bl::Blade)::Multivector
    btm = bitmap(bl)+1
    if typeof(btm) != typeof(bl.algebra.max)
        btm = typeof(bl.algebra.max)(btm)
    end

    values = sparsevec([btm], [scalar(bl)], bl.algebra.max)
    return Multivector(values, bl.algebra)
end

function Multivector(blades::GAVector, al::Algebra = gb_current_algebra)::Multivector
    bitmaps = []
    scalars = []
    for bl in blades
        bl = Blade(bl)
        push!(bitmaps, bitmap(bl)+1)
        push!(scalars, scalar(bl))
    end

    if typeof(al.max) == BigInt
        bitmaps = BigInt.(bitmaps)
    end

    values = sparsevec(bitmaps, scalars, al.max)
    return Multivector(values, al)
end

function str_bl(n::Integer, al::AlgebraFull)::String
    return al.basis_bit_order[n]
end

function str_bl(n::Integer, al::AlgebraMin)::String
    n -= 1
    ans = ""
    #for i in 1:sizeof(n) * 8
    for i in 1:(floor(Int, log2(n))+1)
        if (n >> (i - 1)) & 1 == 1
            ans *= al.symbols[i]
        end
    end
    return ans
end

# Pretty print a GAType
function Base.show(io::IO, mv::GAType)

    if mv.algebra.symbols == []
        print(mv.blade_array)
        return nothing
    end

    bitmaps = mv.blade_array.nzind
    scalars = mv.blade_array.nzval
    ans = ""

    for i in eachindex(bitmaps)
        if bitmaps[i] > 1 && scalars[i] > 0
            ans = ans * "+ " * string(round(scalars[i], digits=3)) * "*" * str_bl(bitmaps[i], mv.algebra) * " "
        elseif bitmaps[i] > 1 && scalars[i] < 0
            ans = ans * "- " * string(round(abs(scalars[i]), digits=3)) * "*" * str_bl(bitmaps[i], mv.algebra) * " "
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

function Base.show(io::IO, ::MIME"text/plain", xs::GAArray)
    print(xs)
end