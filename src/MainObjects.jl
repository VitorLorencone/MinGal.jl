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
end

function Base.length(mv::GAType)::Integer
    dropzeros!(mv.blade_array)
    return length(mv.blade_array.nzind)
end

"""
    Blade(bitmap, scalar)::Blade ||
    Blade(mv::Multivector)::Blade ||
    Blade(k::Number)::Blade

Creates a Blade based on bitmap and Scalars or Converts a multivector into a blade. It may throw an error

# Arguments
- `bitmap::Integer` : The bitmap that specifies what basis vectors are present in this blade 
- `scalar::Number` : The scalar of the basis blade.
||
- `mv::Multivector`
||
- `k::Number`

# Return
Returns a Blade.

"""
function Blade(bitmap::Integer, scalar::Number)
    if typeof(bitmap) != typeof(gb_current_algebra.max)
        bitmap = typeof(gb_current_algebra.max)(bitmap)
    end
    
    values = sparsevec([bitmap+1], [scalar], gb_current_algebra.max)
    return Blade(values)
end

function Blade(mv::Multivector)::Blade
    if length(mv) > 1
        throw(DomainError(mv, "You cannot convert this multivector to a blade."))
    elseif length(mv) == 0
        return Blade(0)
    end
    return Blade(mv.blade_array.nzind[1]-1, mv.blade_array.nzval[1])
end

function Blade(k::Number)::Blade
    bitmap = 1
    
    if typeof(gb_current_algebra.max) == BigInt
        bitmap = typeof(gb_current_algebra.max)(bitmap)
    end
    
    values = sparsevec([bitmap], [k], gb_current_algebra.max)
    return Blade(values)
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
    Multivector(base_vectors, scalars) ||
    Multivector(bl::Blade) ||
    Multivector(blades::Array{Blade})

Constructor function for creating multivectors.

# Arguments
- `base_vectors::Vector{Integer}` : An array of integers, representing the actual basis blade that exists in this object in bit order.
- `scalars::Vector{Number}` : An array of Floats, representing the scalars of each basis blade in bit order.
||
- `bl::Blade` : A Blade, to convert to Multivector.
||
- `blades::GAVector` : An array of Blades, to convert to Multivector.

# Return
Returns a Multivector.

"""
function Multivector(base_vectors::Array, scalars::Array)::Multivector
    if typeof(gb_current_algebra.max) == BigInt
        base_vectors = BigInt.(base_vectors)
    end
    
    values = sparsevec(base_vectors .+ 1, scalars, gb_current_algebra.max)
    return Multivector(values)
end

function Multivector(bl::Blade)::Multivector
    btm = bitmap(bl)+1
    if typeof(btm) != typeof(gb_current_algebra.max)
        btm = typeof(gb_current_algebra.max)(btm)
    end

    values = sparsevec([btm], [scalar(bl)], gb_current_algebra.max)
    return Multivector(values)
end

function Multivector(blades::GAVector)::Multivector
    bitmaps = []
    scalars = []
    for bl in blades
        bl = Blade(bl)
        push!(bitmaps, bitmap(bl)+1)
        push!(scalars, scalar(bl))
    end

    if typeof(gb_current_algebra.max) == BigInt
        bitmaps = BigInt.(bitmaps)
    end

    values = sparsevec(bitmaps, scalars, gb_current_algebra.max)
    return Multivector(values)
end

function str_bl(n::Integer, al::AlgebraFull)::String
    return gb_current_algebra.basis_bit_order[n]
end

function str_bl(n::Integer, al::AlgebraMin)::String
    n -= 1
    ans = ""
    #for i in 1:sizeof(n) * 8
    for i in 1:(floor(Int, log2(n))+1)
        if (n >> (i - 1)) & 1 == 1
            ans *= gb_current_algebra.symbols[i]
        end
    end
    return ans
end

# Pretty print a GAType
function Base.show(io::IO, mv::GAType)

    if gb_current_algebra.symbols == []
        print(mv.blade_array)
        return nothing
    end

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

function Base.show(io::IO, ::MIME"text/plain", xs::GAArray)
    print(xs)
end