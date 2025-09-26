include("GeometricFunctions.jl")

# Functions for operator overloading between types.

# * Geometric Product or Product by Scalar
function Base.:*(mi::GAType, mj::GAType)::GAType
    return geometric_product(mi,mj)
end

function Base.:*(mv::GAType, k::Number)::GAType
    return product_by_scalar(mv,k)
end
function Base.:*(k::Number, mv::GAType)::GAType
    return product_by_scalar(mv,k)
end

# / Division by Scalar
function Base.:/(mv::GAType, k::Number)::GAType
    return product_by_scalar(mv, 1/k)
end

function Base.:/(k::Number, mv::GAType)::Number
    if grade(mv) != 0
        error("You cannot divide by a multivector")
    end
    return k/(mv[0])
end

function Base.:/(mv::GAType, mi::GAType)
    if grade(mv) != 0 && grade(mi) == 0
        return product_by_scalar(mv, 1/(mi[0]))
    elseif grade(mv) == 0 && grade(mi) == 0
        return mv[0]/mi[0]
    else
        error("You cannot divide by a multivector")
    end
end

# \ Inner Product
function Base.:\(mi::GAType, mj::GAType)::GAType
    return inner_product(mi,mj)
end

# ^ Outer Product
function Base.:^(mi::GAType, mj::GAType)::GAType
    return outer_product(mi,mj)
end

# == or != comparison
function Base.:(==)(mi::GAType, mj::GAType)::Bool
    return mi.blade_array == mj.blade_array
end
function Base.:(!=)(mi::GAType, mj::GAType)::Bool
    return mi.blade_array != mj.blade_array
end

# - minus
function Base.:-(mv::GAType)::GAType
    return -1 * mv
end

# + sum
function Base.:+(mi::GAType, mj::GAType)::GAType
    return multivector_sum(mi, mj)
end
function Base.:+(k::Number, mv::GAType)::GAType
    return multivector_sum(Blade(0, k), mv)
end
function Base.:+(mv::GAType, k::Number)::GAType
    return multivector_sum(Blade(0, k), mv)
end

# - sub
function Base.:-(mi::GAType, mj::GAType)::GAType
    return multivector_sub(mi, mj)
end
function Base.:-(k::Number, mv::GAType)::GAType
    return multivector_sub(Blade(0, k), mv)
end
function Base.:-(mv::GAType, k::Number)::GAType
    return multivector_sub(mv, Blade(0, k))
end

# ~ reverse
function Base.:~(mi::GAType)::GAType
    return reverse(mi)
end