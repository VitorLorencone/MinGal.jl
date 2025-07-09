include("BasisFunctions.jl")

"""
    canonical_reordering_sign(ei, ej)::Int

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
    bitmap = bitmap(ei) ⊻ bitmap(ej)
    sign = canonical_reordering_sign(ei, ej)
    scalar = scalar(ei) * scalar(ej) * sign

    meet = bitmap(ei) & bitmap(ej)

    i = 1
    while meet != 0
        if (meet & 1) != 0
            scalar *= gb_current_algebra.metric[i]
        end
        i += 1
        meet = meet >> 1
    end

    return Multivector(Blade(bitmap, scalar))
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
        return Blade(0, 0)
    end

    bitmap = bitmap(ei) ⊻ bitmap(ej)
    sign = canonical_reordering_sign(ei, ej)
    scalar = ei.scalar * ej.scalar * sign

    return Multivector(Blade(bitmap, scalar))
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
    result = mi.blade_array + mj.blade_array
    return Multivector(result.nzind, result.nzval)
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
    result = mi.blade_array - mj.blade_array
    return Multivector(result.nzind, result.nzval)
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
function product_by_scalar(mv::GAType, k::Number)::Number
    return Multivector(mv.nzind, k .* mv.nzval)
end

"""
    geometric_product(ei::GAType, ej::GAType)::GAType

Function that computes the geometric product of two multivectors and return its result.

# Arguments
- `ei::Multivector` : A GAType.
- `ej::Multivector` : A GAType.

# Return
The result GAType.

"""
function geometric_product(ei:: GAType, ej::GAType)::GAType
    result = Multivector([1],[0])
    for (i, si) in zip(ei.blade_array.nzind, ei.blade_array.nzval)
        for (j, sj) in zip(ej.blade_array.nzind, ej.blade_array.nzval)
            result += geometric_product(Blade(i-1, si), Blade(j-1, sj))
        end
    end
    return result
end

"""
    outer_product(ei::GAType, ej::GAType)::GAType

Function that computes the outer product of two multivectors and return its result.

# Arguments
- `ei::Multivector` : A GAType.
- `ej::Multivector` : A GAType.

# Return
The result GAType.

"""
function outer_product(ei:: GAType, ej::GAType)::GAType
    result = Multivector([1],[0])
    for (i, si) in zip(ei.blade_array.nzind, ei.blade_array.nzval)
        for (j, sj) in zip(ej.blade_array.nzind, ej.blade_array.nzval)
            result += outer_product(Blade(i-1, si), Blade(j-1, sj))
        end
    end
    return result
end

"""
    inner_product(ei::GAType, ej::GAType)::GAType

Function that computes the outer product of two multivectors and return its result.

# Arguments
- `ei::Multivector` : A GAType.
- `ej::Multivector` : A GAType.

# Return
The result GAType.

"""
function inner_product(ei:: GAType, ej::GAType)::GAType
    result = Multivector([1],[0])
    for (i, si) in zip(ei.blade_array.nzind, ei.blade_array.nzval)
        for (j, sj) in zip(ej.blade_array.nzind, ej.blade_array.nzval)
            result += inner_product(Blade(i-1, si), Blade(j-1, sj))
        end
    end
    return result
end