include("BasisFunctions.jl")
include("SparseOperationFix.jl")

"""
    canonical_reordering_sign(ei, ej)::Integer

Function that computes the reordering of a Blade to get into canonical order.

# Arguments
- `ei::Blade` : A Blade.
- `ej::Blade` : A Blade.

# Return
The sign value, either 1 or -1.

"""
function canonical_reordering_sign(ei::Blade, ej::Blade)::Int
    i = bitmap(ei)
    j = bitmap(ej)

    i = i >> 1
    sum = 0
    while(i != 0)
        sum = sum + count_ones(i & j)
        i = i >> 1
    end

    return (sum & 1) == 0 ? 1.0 : -1.0
end

"""
    geometric_product(ei, ej)::Multivector

Function that returns the Geometric Product between two blades.

# Arguments
- `ei::Blade` : A Blade.
- `ej::Blade` : A Blade.

# Return
The result Blade.

"""
function geometric_product(ei::Blade, ej::Blade)::Multivector
    bitmap_value = bitmap(ei) ⊻ bitmap(ej)
    sign = canonical_reordering_sign(ei, ej)
    scalar_value = scalar(ei) * scalar(ej) * sign

    meet::Number = bitmap(ei) & bitmap(ej)

    i = 1
    while meet != 0
        if (meet & 1) != 0
            scalar_value *= gb_current_algebra.metric[i]
        end
        i += 1
        meet = meet >> 1
    end

    return Multivector(Blade(bitmap_value, scalar_value))
end

"""
    outer_product(ei, ej)::Multivector

Function that returns the Outer Product between two blades.

# Arguments
- `ei::Blade` : A Blade.
- `ej::Blade` : A Blade.

# Return
The result Blade.

"""
function outer_product(ei::Blade, ej::Blade)::Multivector
    if bitmap(ei) & bitmap(ej) != 0
        return Multivector([0],[0])
    end

    bitmap_value::Number = bitmap(ei) ⊻ bitmap(ej)
    sign::Number = canonical_reordering_sign(ei, ej)
    scalar_value::Number = scalar(ei) * scalar(ej) * sign

    return Multivector(Blade(bitmap_value, scalar_value))
end

"""
    inner_product(ei, ej)::Multivector

Function that returns the Inner Product between two blades.

# Arguments
- `ei::Blade` : A Blade.
- `ej::Blade` : A Blade.

# Return
The result Blade.

"""
function inner_product(ei::Blade, ej::Blade)::Multivector

    k = grade(ei)
    l = grade(ej)
    if k != 0 && l != 0
        return grade_projection(geometric_product(ei, ej), abs(k-l))
    elseif k == 0
        return Multivector(ej)
    else
        return Multivector(ei)
    end
end

"""
    multivector_sum(mi::GAType, mj::GAType)::GAType

Function that sums two GAType and return its result.

# Arguments
- `mi::GAType`
- `mj::GAType``

# Return
The result Multivector.

"""
function multivector_sum(mi::GAType, mj::GAType)::GAType
    result = sparsevec_operate(Base.:+, mi.blade_array, mj.blade_array)
    return Multivector(result)
end

"""
    multivector_sub(mi::GAType, mj::GAType)::GAType

Function that subtracts two multivectors and return its result.

# Arguments
- `mi::GAType`
- `mj::GAType``

# Return
The result GAType.

"""
function multivector_sub(mi::GAType, mj::GAType)::GAType
    result = sparsevec_operate(Base.:-, mi.blade_array, mj.blade_array)
    return Multivector(result)
end

"""
    product_by_scalar(mv::GAType, k::Number)::GAType

Function that calculates the product by scalar between a blade and a scalar.

# Arguments
- `mv::GAType`
- `k::Number``

# Return
The result GAType.

"""
function product_by_scalar(mv::GAType, k::Number)::GAType
    result = sparsevec(mv.blade_array.nzind, k .* mv.blade_array.nzval, gb_current_algebra.max)
    return Multivector(result)
end

"""
    geometric_product(ei::GAType, ej::GAType)::GAType

Function that computes the geometric product of two multivectors and return its result.

# Arguments
- `ei::GAType` : A GAType.
- `ej::GAType` : A GAType.

# Return
The result GAType.

"""
function geometric_product(ei::GAType, ej::GAType)::GAType
    result = Multivector([0],[0])
    for x in ei
        for y in ej
            result += geometric_product(x, y)
        end
    end
    dropzeros!(result.blade_array)
    return result
end

"""
    outer_product(ei::GAType, ej::GAType)::GAType

Function that computes the outer product of two multivectors and return its result.

# Arguments
- `ei::GAType` : A GAType.
- `ej::GAType` : A GAType.

# Return
The result GAType.

"""
function outer_product(ei::GAType, ej::GAType)::GAType    
    result = Multivector([0],[0])
    for x in ei
        for y in ej
            result += outer_product(x, y)
        end
    end
    dropzeros!(result.blade_array)
    return result
end

"""
    inner_product(ei::GAType, ej::GAType)::GAType

Function that computes the outer product of two multivectors and return its result.

# Arguments
- `ei::GAType` : A GAType.
- `ej::GAType` : A GAType.

# Return
The result GAType.

"""
function inner_product(ei::GAType, ej::GAType)::GAType
    result = Multivector([0],[0])
    for x in ei
        for y in ej
            result += inner_product(x, y)
        end
    end
    dropzeros!(result.blade_array)
    return result
end

"""
    reverse(ei::GAType)::GAType

Function that computes the reverse of a multivector and return its result.

# Arguments
- `ei::Multivector` : A GAType.

# Return
The reverse of ei.

"""
function reverse(ei::GAType)::GAType
    result = Multivector([0],[0])
    for x in ei
        k = grade(Blade(i-1, si))
        result += (-1)^(k*(k-1)/2) * x
    end
    return result
end

"""
    involution(ei::GAType)::GAType

Function that computes the involuction of a multivector and return its result.

# Arguments
- `ei::Multivector` : A GAType.

# Return
The involution of ei.

"""
function involution(ei::GAType)::GAType
    result = Multivector([0],[0])
    for x in ei
        k = grade(Blade(i-1, si))
        result += (-1)^(k) * x
    end
    return result
end

"""
    left_contraction(ei::GAType, ej::GAType)::GAType

Function that computes the left contraction of two multivectors and return its result.

# Arguments
- `ei::GAType` : A GAType.
- `ej::GAType` : A GAType.

# Return
The result GAType.

"""
function left_contraction(ei::Blade, ej::Blade)::Multivector
    if grade(ei) <= grade(ej)
        return grade_projection(ei*ej, grade(ej) - grade(ei))
    end
    return Multivector([0],[0])
end

function left_contraction(ei::GAType, ej::GAType)::GAType
    result = Multivector([0],[0])
    for x in ei
        for y in ej
            result += left_contraction(x, y)
        end
    end
    dropzeros!(result.blade_array)
    return result
end

"""
    right_contraction(ei::GAType, ej::GAType)::GAType

Function that computes the right contraction of two multivectors and return its result.

# Arguments
- `ei::GAType` : A GAType.
- `ej::GAType` : A GAType.

# Return
The result GAType.

"""
function right_contraction(ei::Blade, ej::Blade)::Multivector
    if grade(ei) >= grade(ej)
        return grade_projection(ei*ej, grade(ei) - grade(ej))
    end
    return Multivector([0], [0])
end

function right_contraction(ei::GAType, ej::GAType)::GAType
    result = Multivector([0],[0])
    for x in ei
        for y in ej
            result += right_contraction(x, y)
        end
    end
    dropzeros!(result.blade_array)
    return result
end

"""
    grade_selection(ei::GAType, k::Number)::GAType

Function that retrieves the k grade part of a GAType.

# Arguments
- `ei::Multivector` : A GAType.
- `k::Number` : The selected grade.

# Return
The result GAType.

"""
function grade_selection(ei::GAType, k::Number)::GAType
    result = Multivector([0],[0])
    for bl in ei
        if grade(bl) == k
            result += bl
        end
    end
    return result
end

"""
    invert(ei::GAType)::GAType

Function that computes the inverse of a multivector (mostly versors and blades supported!!!!!) and return its result.
May throw an error if not inversible

# Arguments
- `ei::GAType` : A GAType.

# Return
The inverse of ei.

"""
function invert(ei::GAType)::GAType
    ei_rev = ~ei
    denom = ei*ei_rev

    if length(denom) != 1 || denom[0] == 0
        error("Non-invertible GAType")
    end

    return ei_rev/denom[0]
end

"""
    dual(ei::GAType)::GAType

Function that computes the dual of a multivector and return its result.
Just be careful with r>1.

# Arguments
- `ei::Multivector` : A GAType.

# Return
The dual of ei.

"""
function dual(ei::GAType)::GAType
    if gb_current_algebra.r == 0
        key = gb_current_algebra.max-1
        I = Multivector([key], [1])
        return ei*I*canonical_reordering_sign(Blade(key, 1), Blade(key, 1)) 
    else
        result = Multivector([0], [0])
        n = gb_current_algebra.max
        for i in ei
            key = bitmap(i)
            value = scalar(i)
            dual_index = n - key - 1
            sign = canonical_reordering_sign(Blade(key, 1), Blade(dual_index, 1))
            if sign < 0
                result[dual_index] = -value
            else
                result[dual_index] = value
            end
        end
        return result
    end
end

"""
    dual(ei::GAType)::GAType

Function that computes the undual of a multivector and return its result.
Just be careful with r>1.

# Arguments
- `ei::Multivector` : A GAType.

# Return
The undual of ei.

"""
function undual(ei::GAType)::GAType
    I = Multivector([gb_current_algebra.max-1], [1])

    if gb_current_algebra.r == 0
        return ei<<I

    else
        result = Multivector([0], [0])
        n = gb_current_algebra.max
        for i in ei
            key = bitmap(i)
            value = scalar(i)
            dual_index = n - key - 1
            sign = canonical_reordering_sign(Blade(dual_index, 1), Blade(key, 1))
            if sign < 0
                result[dual_index] = -value
            else
                result[dual_index] = value
            end
        end
        return result
    end
end

"""
    regressive_product(ei::GAType, ej::GAType)::GAType

Function that computes the regressive product of two multivectors and return its result.

# Arguments
- `ei::GAType` : A GAType.
- `ej::GAType` : A GAType.

# Return
The result GAType.

"""
function regressive_product(ei::Blade, ej::Blade)::Multivector
    return undual(dual(ei)^dual(ej))
end

function regressive_product(ei::GAType, ej::GAType)::GAType
    result = Multivector([0],[0])
    for x in ei
        for y in ej
            result += regressive_product(Blade(i-1, si), Blade(j-1, sj))
        end
    end
    dropzeros!(result.blade_array)
    return result
end

"""
    isScalar(ei::GAType)::Bool

Function that checks where a GAType is a Scalar or not.

# Arguments
- `ei::GAType` : A GAType.

"""
function isScalar(ei::GAType)::Bool
    dropzeros!(ei.blade_array)
    return ei.blade_array.nzind == [1] || length(ei) == 0
end

"""
    isBlade(ei::GAType)::Bool

Function that checks where a GAType is a Blade or not.

# Arguments
- `ei::GAType` : A GAType.

"""
function isBlade(ei::GAType)::Bool
    dropzeros!(ei.blade_array)
    return length(ei) == 1
end

# TODO generic exp whe A² is not a scalar
"""
    exp_ga(ei::GAType)::GAType

Function that computes the exponential of elements that squares to a scalar.
A more generic function was not yet implemented

# Arguments
- `ei::GAType` : A GAType.

# Return
The result GAType.

"""
function exp_ga(ei::GAType)::GAType
    A2 = ei*ei
    if isScalar(A2) && A2[0] == 0
        return ei+1
    elseif isScalar(A2)
        k = A2[0]
        if k < 0
            a = sqrt(-k)
            return (sin(a)/a) * ei + cos(a)
        else
            a = sqrt(k)
            return (sinh(a)/a) * ei + cosh(a)
        end
    else
        error("There is no implementation of exp when element do not squares to a scalar")
    end
end