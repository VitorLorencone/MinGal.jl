include("BasisFunctions.jl")
include("SparseOperationFix.jl")

"""
    canonical_reordering_sign(ei, ej)::Integer

Function that computes the reordering of a Blade to get into canonical order.

# Arguments
- `ei::Blade`
- `ej::Blade`

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
- `ei::Blade`
- `ej::Blade`

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
- `ei::Blade`
- `ej::Blade`

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
- `ei::Blade`
- `ej::Blade`

# Return
The result Blade.

"""
function inner_product(ei::Blade, ej::Blade)::Multivector
    k = grade(ei)
    l = grade(ej)

    if k != 0 && l != 0
        return grade_projection(geometric_product(ei, ej), abs(k-l))
    else
        return geometric_product(ei, ej)
    end
end

"""
    multivector_sum(mi::GAType, mj::GAType)::GAType

Function that sums two GAType and return its result.

# Arguments
- `mi::GAType`
- `mj::GAType`

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
- `mj::GAType`

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
- `k::Number`

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
- `ei::GAType`
- `ej::GAType`

# Return
The result GAType.

"""
function geometric_product(ei::GAType, ej::GAType)::GAType
    result = Multivector([0],[0])
    # x and y are both Blade types
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
- `ei::GAType`
- `ej::GAType`

# Return
The result GAType.

"""
function outer_product(ei::GAType, ej::GAType)::GAType    
    result = Multivector([0],[0])
    # x and y are both Blade types
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
- `ei::GAType`
- `ej::GAType`

# Return
The result GAType.

"""
function inner_product(ei::GAType, ej::GAType)::GAType
    result = Multivector([0],[0])
    # x and y are both Blade types
    for x in ei
        for y in ej
            result += inner_product(x, y)
        end
    end
    dropzeros!(result.blade_array)
    return result
end

"""
    revert(ei::GAType)::GAType

Function that computes the reverse of a multivector and return its result.

# Arguments
- `ei::Multivector`

# Return
The reverse of ei.

"""
function revert(ei::GAType)::GAType
    result = Multivector([0],[0])
    for x in ei
        k = grade(x)
        result += (-1)^(k*(k-1)/2) * x
    end
    return result
end

"""
    conjugate(ei::GAType)::GAType

Function that computes the conjugate of a multivector and return its result.

# Arguments
- `ei::Multivector`

# Return
The conjugate of ei.

"""
function conjugate(ei::GAType)::GAType
    result = Multivector([0],[0])
    for x in ei
        k = grade(x)
        r = grade_minus(x)
        result += (-1)^r * (-1)^(k*(k-1)/2) * x
    end
    return result
end

"""
    clifford_conjugation(ei::GAType)::GAType

Function that computes the Clifford Conjugation of a multivector and return its result.

# Arguments
- `ei::Multivector`

# Return
The Clifford Conjugate of ei.

"""
function clifford_conjugation(ei::GAType)::GAType
    result = Multivector([0],[0])
    for x in ei
        k = grade(x)
        result += (-1)^(k*(k+1)/2) * x
    end
    return result
end

"""
    grade_involution(ei::GAType)::GAType

Function that computes the grade involuction of a multivector and return its result.

# Arguments
- `ei::Multivector`

# Return
The grade involution of ei.

"""
function grade_involution(ei::GAType)::GAType
    result = Multivector([0],[0])
    for x in ei
        k = grade(x)
        result += (-1)^(k) * x
    end
    return result
end

"""
    left_contraction(ei::GAType, ej::GAType)::GAType

Function that computes the left contraction of two multivectors and return its result.

# Arguments
- `ei::GAType`
- `ej::GAType`

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
- `ei::GAType`
- `ej::GAType`

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

Function that retrieves the k grade part of a GAType. Duplicated, since grade_projection 
does the same, but a new algorithm and more generalized.

# Arguments
- `ei::Multivector`
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

Function that computes the inverse of a multivector (mostly versors and blades supported!) and return its result.
May throw an error if not invertible.

# Arguments
- `ei::GAType`

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
Just be careful with parameter r > 0. It uses Polarity algorithm for
invertible algebras and hodge dual for degenerate algebras.

# Arguments
- `ei::Multivector`

# Return
The dual of ei.

"""
function dual(ei::GAType)::GAType
    if gb_current_algebra.r == 0
        # Polarity
        key = gb_current_algebra.max-1
        I = Multivector([key], [1])
        return ei*I*canonical_reordering_sign(Blade(key, 1), Blade(key, 1)) 
    else
        # Hodge
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
    undual(ei::GAType)::GAType

Function that computes the undual of a multivector and return its result.
Just be careful with r > 0. It uses Polarity algorithm for invertible 
algebras and hodge dual for degenerate algebras.

# Arguments
- `ei::Multivector`

# Return
The undual of ei.

"""
function undual(ei::GAType)::GAType
    if gb_current_algebra.r == 0
        # Polarity
        I = Multivector([gb_current_algebra.max-1], [1])
        return ei<<I
    else
        # Hodge
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
- `ei::GAType`
- `ej::GAType`

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
            result += regressive_product(x, y)
        end
    end
    dropzeros!(result.blade_array)
    return result
end

"""
    isScalar(ei::GAType)::Bool

Function that checks where a GAType is a Scalar or not.

# Arguments
- `ei::GAType`

# Return
True if it's a Scalar, False otherwise.

"""
function isScalar(ei::GAType)::Bool
    dropzeros!(ei.blade_array)
    return ei.blade_array.nzind == [1] || length(ei) == 0
end

"""
    isBlade(ei::GAType)::Bool

Function that checks where a GAType is a Blade or not.

# Arguments
- `ei::GAType`

# Return
True if it's a Scalar, False otherwise.

"""
function isBlade(ei::GAType)::Bool
    dropzeros!(ei.blade_array)
    return length(ei) == 1
end

# TODO generic exp when A² is not a scalar
"""
    exp_ga(ei::GAType)::GAType

Function that computes the exponential of elements that squares to a scalar.
A more generic function was not yet implemented.

# Arguments
- `ei::GAType`

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

"""
    norm(ei::GAType)::Number

Function that computes the reverse norm in absolute value (conjugate norm) for a GAType on usual algebras. 
For degenerate algebras, computes euclidian norm.

# Arguments
- `ei::GAType`

# Return
The norm of ei.

"""
function norm(ei::GAType)::Number
    if gb_current_algebra.r == 0
        return conjugate_norm(ei)
    else
        return euclidian_norm(ei)
    end
end

"""
    reverse_norm(ei::GAType)::Number

Function that computes the reverse norm of a GAType on usual algebras. In practice 
the reverse norm is useful, especially due to its possible negative sign

# Arguments
- `ei::GAType`

# Return
The reverse norm of ei.

"""
function reverse_norm(ei::GAType)::Number
    rev_norm = (ei*revert(ei))[0]
    return sign(rev_norm) * sqrt(abs(rev_norm))
end

"""
    euclidian_norm(ei::GAType)::Number

Function that computes the euclidian coefficients norm of a GAType. It is also good
for degenerate algebras then the metric might be ignored.

# Arguments
- `ei::GAType`

# Return
The euclidian norm of ei.

"""
function euclidian_norm(ei::GAType)::Number
    sum = 0
    for x in ei
        sum += scalar(x) * scalar(x)
    end
    return sqrt(sum)
end

"""
    conjugate_norm(ei::GAType)::Number

Function that computes the conjugate norm of a GAType. Also another possible norm
for degenerate algebras.

# Arguments
- `ei::GAType`

# Return
The conjugate norm of ei.

"""
function conjugate_norm(ei::GAType)::Number
    sqrt(abs((ei*conjugate(ei))[0]))
end

"""
    clifford_norm(ei::GAType)::Number

Function that computes the clifford norm of a GAType.

# Arguments
- `ei::GAType`

# Return
The clifford norm of ei.

"""
function clifford_norm(ei::GAType)::Number
    sqrt(abs((ei*clifford_conjugation(ei))[0]))
end