include("MainObjects.jl")

"""
    grade(bl::Blade)::Int
    grade(mv::Multivector)::Int

Function that returns the grade of the Blade or multivector.

# Arguments
- `bl::Blade`
- `mv::Multivector`

# Return
An integer, the grade of the blade or multivector.

"""
function grade(bl::Blade)::Int
    return count_ones(bitmap(bl))
end

function grade(mv::Multivector)::Int
    grade::Int = 0
    for bitmap in mv.nzind
        grade = max(grade, bitmap-1)
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
- `k::Int` : An integer to the Grade Projection
- `mv::Blade` : A Multivector.

# Return
The result Blade. It might be the 1D blade "1"

"""
function grade_projection(bl::Blade, k::Int)::Blade
    if grade(bl) == k
        return Multivector(bl)
    else
        return Multivector(Blade(0,1))
    end
end

function grade_projection(mv::Multivector, k::Int)::Blade
    bl = Blade(mv)

    if grade(bl) == k
        return Multivector(bl)
    else
        return Multivector(Blade(0,1))
    end
end

"""
    scalar_product(ei, ej)::Int

Function that returns the Scalar Product between two basis blades.

# Arguments
- `ei::Blade` : A Blade.
- `ej::Blade` : A Blade.

# Return
The result Integer.

"""
function scalar_product(ei::GAType, ej::GAType)::Int
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
- `k::Int` : The scalar.
||
- `ei::Blade` : A Blade.

# Return
The scalar value.

"""
function get_scalar(mv::Multivector, k::Int)::Number
    if k in mv.blade_array.nzind
        return mv.blade_array[k]
    end
    return 0
end

function get_scalar(mv::Multivector, ei::Blade)::Number
    if (bitmap(ei)+1) in mv.blade_array.nzind
        return mv.blade_array[bitmap(ei)+1]
    end
    return 0
end