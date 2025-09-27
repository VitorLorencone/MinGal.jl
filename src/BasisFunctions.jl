include("MainObjects.jl")

"""
    grade(bl::Blade)::Integer
    grade(mv::Multivector)::Integer

Function that returns the grade of the Blade or multivector.

# Arguments
- `bl::Blade`
- `mv::Multivector`

# Return
An integer, the grade of the blade or multivector.

"""
function grade(bl::Blade)::Integer
    return count_ones(bitmap(bl))
end

function grade(mv::Multivector)::Integer
    grade::Integer = 0
    for bitmap in mv.blade_array.nzind
        grade = max(grade, count_ones(bitmap-1))
    end
    return grade
end

"""
    grade_projection(bl, k)::Blade ||
    grade_projection(mv, k)::Blade

Function that returns the grade Projection between a Blade and an Integer. In the
multivector form, it requires that it has only one blade.

# Arguments
- `bl::Blade` : A Blade.
- `k::Integer` : An integer to the Grade Projection
- `mv::Blade` : A Multivector.

# Return
The result Blade. It might be the 1D blade "1"

"""
function grade_projection(bl::Blade, k::Integer)::Multivector
    if grade(bl) == k
        return Multivector(bl)
    else
        return Multivector(Blade(0,0))
    end
end

function grade_projection(mv::Multivector, k::Integer)::Multivector
    bl = Blade(mv)

    if grade(bl) == k
        return Multivector(bl)
    else
        return Multivector([0],[0])
    end
end

"""
    scalar_product(ei, ej)::Integer

Function that returns the Scalar Product between two basis blades.

# Arguments
- `ei::Blade` : A Blade.
- `ej::Blade` : A Blade.

# Return
The result Integer.

"""
function scalar_product(ei::GAType, ej::GAType)::Integer
    ei = Blade(ei) # Checks if it is a Blade
    ej = Blade(ej) # Checks if it is a Blade

    if(grade(ei) != 1 || scalar(ei) == 0)
        throw(DomainError(ei, "This operation must be executed with basis blades."))
    elseif(grade(ej) != 1 || scalar(ej) == 0)
        throw(DomainError(ej, "This operation must be executed with basis blades."))
    end

    i = log2(bitmap(ei))
    j = log2(bitmap(ej))
    
    if i != j
        return 0
    elseif i >= 1 && j <= gb_current_algebra.r
        return 0
    elseif i > gb_current_algebra.r && j <= gb_current_algebra.r + gb_current_algebra.p
        return 1
    else
        return -1
    end
end

"""
    get_scalar(mv, k)::Number ||
    get_scalar(mv, ei)::Number

Function that returns the scalar value in index k from multivector mv.
The index k follows the bit order basis.The ei value represents the 
blade that you want the scalar from.

# Arguments
- `mv::Multivector` : A Multivector.
- `k::Integer` : The scalar.
||
- `ei::Blade` : A Blade.

# Return
The scalar value.

"""
function get_scalar(mv::Multivector, k::Integer)::Number
    if (k+1) in mv.blade_array.nzind
        return mv.blade_array[k+1]
    end
    return 0
end

function get_scalar(mv::Multivector, ei::Multivector)::Number
    bl = Blade(ei)
    if (bitmap(bl)+1) in mv.blade_array.nzind
        return mv.blade_array[bitmap(bl)+1]
    end
    return 0
end

"""
    set_scalar(mv, val, k)::Number ||
    set_scalar(mv, val, ei)::Number

Function that changes the scalar value in index k from multivector mv.
The index k follows the bit order basis.The ei value represents the 
blade that you want to change the scalar from.

# Arguments
- `mv::Multivector`
- `val::Number` : New Value
- `k::Integer` : Index
||
- `ei::Blade` : Index

# Return
The scalar value.

"""
function set_scalar(mv::GAType, val::Number, k::Integer)
    mv.blade_array[k+1] = val
    return val
end

function set_scalar(mv::GAType, val::Number, ei::Multivector)
    bl = Blade(ei)
    mv.blade_array[bitmap(bl)+1] = val
    return val
end

"""
    has_key(mv, k)::Bool ||
    has_key(mv, ei)::Bool

Function that checks if the multivector has a blade.

# Arguments
- `mv::Multivector`
- `k::Integer` : Index
||
- `ei::Blade` : Index

# Return
True, if yes, False, if not.

"""
function has_key(mv::GAType, k::Integer)::Bool
    if(get_scalar(mv, k) != 0)
        return true
    end
    return false
end

function has_key(mv::GAType, ei::Multivector)::Bool
    if(get_scalar(mv, ei) != 0)
        return true
    end
    return false
end

# Base functions Substitution

function Base.getindex(mv::GAType, k::Integer)
    return get_scalar(mv, k)
end

function Base.getindex(mv::GAType, ei::Multivector)
    return get_scalar(mv, ei)
end

function Base.setindex!(mv::GAType, val::Number, k::Integer)
    return set_scalar(mv, val, k)
end

function Base.setindex!(mv::GAType, val::Number, ei::Multivector)
    return set_scalar(mv, val, ei)
end

function Base.iterate(mv::GAType, state = 1)
    inds = mv.blade_array.nzind
    vals = mv.blade_array.nzval

    if state > length(inds)
        return nothing
    else
        return (Blade(inds[state]-1, vals[state]), state + 1)
    end
end

function Base.size(mv::GAType)
    return mv.blade_array.n
end

function Base.keys(mv::GAType)
    return mv.blade_array.nzind .- 1
end

function Base.values(mv::GAType)
    return mv.blade_array.nzval
end

function Base.haskey(mv::GAType, k::Integer)
    return has_key(mv, k)
end

function Base.haskey(mv::GAType, ei::Multivector)
    return has_key(mv, ei)
end