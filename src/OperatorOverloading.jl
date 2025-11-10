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
function Base.:\(mi::GAType, mj::Number)::GAType
    return inner_product(mi, Blade(mj))
end
function Base.:\(mi::Number, mj::GAType)::GAType
    return inner_product(Blade(mi), mj)
end

# ^ Outer Product
function Base.:^(mi::GAType, mj::GAType)::GAType
    return outer_product(mi,mj)
end
function Base.:^(mi::Number, mj::GAType)::GAType
    return inner_product(Blade(mi), mj)
end

# ^ Exponentiation
function Base.:^(mi::GAType, mj::Number)::GAType
    result = Multivector([0],[1])
    for _ in 1:mj
        result = geometric_product(result, mi)
    end
    return result
end

# ^-n invert function
function Base.inv(mi::GAType)::GAType
    return invert(mi)
end

# == or != comparison
function Base.:(==)(mi::GAType, mj::GAType)::Bool
    return mi.blade_array == mj.blade_array
end
function Base.:(==)(mi::GAType, mj::Number)::Bool
    return mi.blade_array == Blade(mj).blade_array
end
function Base.:(==)(mi::Number, mj::GAType)::Bool
    return Blade(mi).blade_array == mj.blade_array
end

function Base.:(!=)(mi::GAType, mj::GAType)::Bool
    return mi.blade_array != mj.blade_array
end
function Base.:(!=)(mi::GAType, mj::Number)::Bool
    return mi.blade_array != Blade(mj).blade_array
end
function Base.:(!=)(mi::Number, mj::GAType)::Bool
    return Blade(mi).blade_array != mj.blade_array
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
    return revert(mi)
end

# << left contraction
function Base.:(<<)(mi::GAType, mj::GAType)::GAType
    left_contraction(mi, mj)
end
function Base.:(<<)(mi::Number, mj::GAType)::GAType
    left_contraction(Blade(mi), mj)
end
function Base.:(<<)(mi::GAType, mj::Number)::GAType
    left_contraction(mi, Blade(mj))
end

# >> right contraction
function Base.:(>>)(mi::GAType, mj::GAType)::GAType
    right_contraction(mi, mj)
end
function Base.:(>>)(mi::Number, mj::GAType)::GAType
    right_contraction(Blade(mi), mj)
end
function Base.:(>>)(mi::GAType, mj::Number)::GAType
    right_contraction(mi, Blade(mj))
end

# & regressive product
function Base.:&(mi::GAType, mj::GAType)::GAType
    regressive_product(mi, mj)
end
function Base.:&(mi::Number, mj::GAType)::GAType
    regressive_product(Blade(mi), mj)
end
function Base.:&(mi::GAType, mj::Number)::GAType
    regressive_product(mi, Blade(mj))
end