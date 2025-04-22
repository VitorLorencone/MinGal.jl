include("OperationTable.jl")

# Functions for operator overloading between types.

function Base.:*(ei::Blade, ej::Blade)::Blade
    return bladeGeometricProduct(ei,ej)
end

function Base.:*(ei::Blade, k::Number)::Blade
    return bladeScalarProduct(ei,k)
end
function Base.:*(k::Number, ei::Blade)::Blade
    return bladeScalarProduct(ei,k)
end

function Base.:*(k::Number, ei::Multivector)::Multivector
    return multivectorByScalar(ei, k)
end
function Base.:*(ei::Multivector, k::Number)::Multivector
    return multivectorByScalar(ei, k)
end
function Base.:/(ei::Multivector, k::Number)::Multivector
    return multivectorByScalar(ei, 1/k)
end

function Base.:\(ei::Blade, ej::Blade)::Blade
    return bladeInnerProduct(ei,ej)
end

function Base.:^(ei::Blade, ej::Blade)::Blade
    return bladeOuterProduct(ei,ej)
end

function Base.:(==)(ei::AbstractGeometricAlgebraType, ej::AbstractGeometricAlgebraType)::Bool
    return ei.val == ej.val
end
function Base.:(!=)(ei::AbstractGeometricAlgebraType, ej::AbstractGeometricAlgebraType)::Bool
    return !(ei == ej)
end

function Base.:-(ei::AbstractGeometricAlgebraType)::AbstractGeometricAlgebraType
    return -1 * ei
end

function Base.:+(ei::AbstractGeometricAlgebraType, ej::AbstractGeometricAlgebraType)::AbstractGeometricAlgebraType
    return multivectorSum(ei, ej)
end
function Base.:+(ei::Int, ej::AbstractGeometricAlgebraType)::AbstractGeometricAlgebraType
    return multivectorSum(Multivectors([1],[ei]), ej)
end
function Base.:+(ei::AbstractGeometricAlgebraType, ej::Int)::AbstractGeometricAlgebraType
    return multivectorSum(ei, Multivectors([1],[ej]))
end

function Base.:-(ei::AbstractGeometricAlgebraType, ej::AbstractGeometricAlgebraType)::AbstractGeometricAlgebraType
    return multivectorSub(ei, ej)
end
function Base.:-(ei::Int, ej::AbstractGeometricAlgebraType)::AbstractGeometricAlgebraType
    return multivectorSub(Multivectors([1],[ei]), ej)
end
function Base.:-(ei::AbstractGeometricAlgebraType, ej::Int)::AbstractGeometricAlgebraType
    return multivectorSub(ei, Multivectors([1],[ej]))
end

function Base.:*(ei::Multivector, ej::Multivector)::AbstractGeometricAlgebraType
    return multivectorGP(ei, ej)
end
function Base.:*(ei::Blade, ej::Multivector)::AbstractGeometricAlgebraType
    return multivectorGP(ei, ej)
end
function Base.:*(ei::Multivector, ej::Blade)::AbstractGeometricAlgebraType
    return multivectorGP(ei, ej)
end

function Base.:\(ei::Multivector, ej::Multivector)::AbstractGeometricAlgebraType
    return multivectorIP(ei, ej)
end
function Base.:\(ei::Blade, ej::Multivector)::AbstractGeometricAlgebraType
    return multivectorIP(ei, ej)
end
function Base.:\(ei::Multivector, ej::Blade)::AbstractGeometricAlgebraType
    return multivectorIP(ei, ej)
end

function Base.:^(ei::Multivector, ej::Multivector)::AbstractGeometricAlgebraType
    return multivectorOP(ei, ej)
end
function Base.:^(ei::Blade, ej::Multivector)::AbstractGeometricAlgebraType
    return multivectorOP(ei, ej)
end
function Base.:^(ei::Multivector, ej::Blade)::AbstractGeometricAlgebraType
    return multivectorOP(ei, ej)
end