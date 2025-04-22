include("MainStructs.jl")

"""
    bladeindex(vec, Al)::Array

Function that returns the indexes of a blade.

# Arguments
- `vec::Blade` : A Blade.
- `Al::AlgebraStruct` : The Algebra, it is setted as CurrentAlgebra.

# Return
Returns an array with all indexes of that blade.

"""
function bladeindex(vec::Blade, Al::AlgebraStruct = CurrentAlgebra)::Array{}
    return Al.Indexes[vec.val.nzind]
end

"""
    bladescalar(vec)::Any

Function that returns scalar of a Blade.

# Arguments
- `vec::Blade` : A Blade.

# Return
A real number, the scalar.

"""
function bladescalar(vec::Blade)::Any
    return vec.val.nzval[1]
end

"""
    lenElements(vec)::Int

Function that returns the length of elements in an Abstract Geometric Algebra Type.

# Arguments
- `vec::AbstractGeometricAlgebraType` : A multivector or a Blade.

# Return
An integer, the ammount of elements.

"""
function lenElements(vec::AbstractGeometricAlgebraType)::Int
    _count_ = 0
    for i in vec.val.nzval
        if(i != 0)
            _count_ += 1
        end
    end
    return _count_
end

"""
    grade(vec)::Int

Function that returns the grade of the Blade.

# Arguments
- `vec::Blade` : A Blade.

# Return
An integer, the grade of the blade.

"""
function grade(vec::Blade)::Int
    if bladeindex(vec)[1] == [0]
        return 0
    end
    return length(bladeindex(vec)[1])
end

"""
    gradeprojection(vec, k)::Blade

Function that returns the grade Projection between a Blade and an Integer.

# Arguments
- `vec::Blade` : A Blade.
- `k::Int` : An integer to the Grade Projection

# Return
The result Blade. It might be the 1D blade "1"

"""
function gradeprojection(vec::Blade, k::Int)::Blade
    if grade(vec) == k
        return vec
    else
        return Multivectors([1],[0])
    end
end

"""
    scalarproduct(ei, ej, Al)::Int

Function that returns the Scalar Product between two basis blades.

# Arguments
- `ei::Blade` : A Blade.
- `ej::Blade` : A Blade.
- `Al::AlgebraStruct` : The Algebra, it is setted as CurrentAlgebra.

# Return
The result Integer.

"""
function scalarproduct(ei::Blade, ej::Blade, Al::AlgebraStruct = CurrentAlgebra)::Int
    i = bladeindex(ei)[1]
    j = bladeindex(ej)[1]

    if(length(i) != 1 || length(j) != 1 || i[1] == 0 || j[1] == 0)
        throw(DomainError(i, "This operation must be executed with basis blades."))
    end
    
    if i[1] != j[1]
        return 0
    elseif i[1] >= 1 && j[1] <= Al.p
        return 1
    else
        return -1
    end
end

"""
    bladeScalarProduct(ei, k)::Blade

Function that returns the Scalar Product between a blades and a 1D Number.

# Arguments
- `ei::Blade` : A Blade.
- `k::Number` : A scalar.

# Return
The result Blade.

"""
function bladeScalarProduct(ei::Blade, k::Number)::Blade
    finalScalar = ei.val.nzval[1] * k
    return Multivectors([ei.val.nzind[1]], [finalScalar])
end

"""
    getScale(ei, k)::Number

Function that returns the Scale value in index k from multivector ei.
The index k follows the ordered set of basis, the same shown in AlgebraStruct.Basis

# Arguments
- `ei::Multivector` : A Multivector.
- `k::Number` : A scalar.

# Return
The Scale value.

"""
function getScale(ei::Multivector, k::Number)::Number
    return ei.val.nzval[k]
end


"""
    getScale(ei, symbol)::Number

Function that returns the Scale value in index k from multivector ei.
The symbol value represents the blade that you want the scale from

# Arguments
- `ei::Multivector` : A Multivector.
- `symbol::String` : A String.

# Return
The Scale value.

"""
function getScale(ei::Multivector, symbol::String)::Number
    index = findfirst(==(symbol), [tupla[1] for tupla in Mingal.CurrentAlgebra.Basis])
    return ei.val.nzval[index]
end