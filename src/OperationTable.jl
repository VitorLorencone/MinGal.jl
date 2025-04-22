include("GaFunctions.jl")

"""
    rawBladeGeometricProduct(ei, ej, Al)::Blade

Function that returns the Geometric Product between two blades.
It is used for the Operation Table, it is high cost for single operations over time.

# Arguments
- `ei::Blade` : A Blade.
- `ej::Blade` : A Blade.
- `Al::Algebra` : The Algebra, it is setted as CurrentAlgebra.

# Return
The result Blade.

"""
function rawBladeGeometricProduct(ei::Blade, ej::Blade, Al::Algebra = CurrentAlgebra)::Blade
    i = bladeindex(ei)
    j = bladeindex(ej)

    finalScalar = bladescalar(ei) * bladescalar(ej)
    finalIndex = setdiff(union(i[1], j[1]), intersect(i[1], j[1]))

    if i[1][1] == 0 && j[1][1] == 0
        return Multivectors([1],[finalScalar])
    elseif i[1][1] == 0
        return bladeScalarProduct(ej, bladescalar(ei))
    elseif j[1][1] == 0
        return bladeScalarProduct(ei, bladescalar(ej))
    end
        
    if i[1][1] == 0
        return ej
    elseif j[1][1] == 0
        return ei
    end

    function signalSwitch(i, j):Int
        if length(i[1]) == 0 || length(j[1]) == 0
            return 1
        end
        for k in eachindex(i[1])
            for l in eachindex(j[1])
                if i[1][k] == j[1][l]
                    if i[1][k] >= 1 && j[1][l] <= Al.p
                        resp = (-1)^(length(i[1]) - k + l - 1)
                        i[1] = filter(x -> x != i[1][k], i[1])
                        j[1] = filter(x -> x != j[1][l], j[1])
                        return resp * signalSwitch(i, j)
                    else
                        resp = (-1)^(length(i[1]) - k + l)
                        i[1] = filter(x -> x != i[1][k], i[1])
                        j[1] = filter(x -> x != j[1][l], j[1])
                        return resp * signalSwitch(i, j)
                    end
                end
            end
        end
        return 1
    end

    finalScalar = finalScalar*signalSwitch(i, j)
    
    if length(finalIndex) == 0
        finalIndex = [0]
    end

    swap = 0
    for i in 2:length(finalIndex)
        key = finalIndex[i]
        j = i - 1
        while j > 0 && finalIndex[j] > key
            finalIndex[j + 1] = finalIndex[j]
            swap+=1
            j -= 1
        end
        finalIndex[j + 1] = key
    end

    finalScalar *= (-1)^(swap)

    return Multivectors([findfirst(x -> x == finalIndex, Al.Indexes)],[finalScalar])

end

"""
    CreateGPTable(Al::Algebra)

Function that creates the Operation Table for the Geometric Product.

# Arguments
- `Al::Algebra` : The Algebra, it is setted as CurrentAlgebra.

"""
function CreateGPTable(Al::Algebra = CurrentAlgebra)
    size = 2^(Al.p + Al.q)
    Table = Array{Any,2}(undef,size,size)

    for i in 1:size
        for j in 1:size
            Table[i,j] = rawBladeGeometricProduct(Multivectors([i],[1]), Multivectors([j],[1]))
        end
    end

    global GPTable = Table
end

"""
    rawBladeInnerProduct(ei, ej)::Blade

Function that returns the Inner Product between two blades.
It is used for the Operation Table, it is high cost for single operations over time.

# Arguments
- `ei::Blade` : A Blade.
- `ej::Blade` : A Blade.

# Return
The result Blade.

"""
function rawBladeInnerProduct(ei::Blade, ej::Blade, GPTable::Array{Any,2} = CreateGPTable())::Blade

    k = grade(ei)
    l = grade(ej)
    if k != 0 && l != 0
        return gradeprojection(GPTable[ei.val.nzind[1], ej.val.nzind[1]], abs(k-l))
    elseif k == 0
        return ej
    else
        return ei
    end

end

"""
    rawBladeOuterProduct(ei, ej)::Blade

Function that returns the Outer Product between two blades.
It is used for the Operation Table, it is high cost for single operations over time.

# Arguments
- `ei::Blade` : A Blade.
- `ej::Blade` : A Blade.

# Return
The result Blade.

"""
function rawBladeOuterProduct(ei::Blade, ej::Blade, GPTable::Array{Any,2} = CreateGPTable())::Blade
    k = grade(ei)
    l = grade(ej)
    return gradeprojection(GPTable[ei.val.nzind[1], ej.val.nzind[1]], l+k)
end

"""
    CreateIPTable(Al::Algebra)

Function that creates the Operation Table for the Inner Product.

# Arguments
- `Al::Algebra` : The Algebra, it is setted as CurrentAlgebra.

"""
function CreateIPTable(GPTable::Array{Any,2} = GPTable, Al::Algebra = CurrentAlgebra)
    size = 2^(Al.p + Al.q)
    Table = Array{Any,2}(undef,size,size)

    for i in 1:size
        for j in 1:size
            Table[i,j] = rawBladeInnerProduct(Multivectors([i],[1]), Multivectors([j],[1]), GPTable)
        end
    end

    global IPTable = Table
end

"""
    CreateOPTable(Al::Algebra)

Function that creates the Operation Table for the Outer Product.

# Arguments
- `Al::Algebra` : The Algebra, it is setted as CurrentAlgebra.

"""
function CreateOPTable(GPTable::Array{Any,2} = GPTable, Al::Algebra = CurrentAlgebra)
    size = 2^(Al.p + Al.q)
    Table = Array{Any,2}(undef,size,size)

    for i in 1:size
        for j in 1:size
            Table[i,j] = rawBladeOuterProduct(Multivectors([i],[1]), Multivectors([j],[1]), GPTable)
        end
    end

    global OPTable = Table
end

"""
    bladeGeometricProduct(ei, ej, Al)::Blade

Function that returns the Geometric Product between two blades with the Operation Table.

# Arguments
- `ei::Blade` : A Blade.
- `ej::Blade` : A Blade.
- `Al::Algebra` : The Algebra, it is setted as CurrentAlgebra.

# Return
The result Blade.

"""
function bladeGeometricProduct(ei::Blade, ej::Blade, Al::Algebra = CurrentAlgebra)::Blade
    return GPTable[ei.val.nzind[1], ej.val.nzind[1]] * bladescalar(ei) * bladescalar(ej)
end

"""
    bladeInnerProduct(ei, ej)::Blade

Function that returns the Inner Product between two blades with the Operation Table.

# Arguments
- `ei::Blade` : A Blade.
- `ej::Blade` : A Blade.

# Return
The result Blade.

"""
function bladeInnerProduct(ei::Blade, ej::Blade)::Blade
    return IPTable[ei.val.nzind[1], ej.val.nzind[1]] * bladescalar(ei) * bladescalar(ej)
end

"""
    bladeOuterProduct(ei, ej)::Blade

Function that returns the Outer Product between two blades with the Operation Table.

# Arguments
- `ei::Blade` : A Blade.
- `ej::Blade` : A Blade.

# Return
The result Blade.

"""
function bladeOuterProduct(ei::Blade, ej::Blade)::Blade
    return OPTable[ei.val.nzind[1], ej.val.nzind[1]] * bladescalar(ei) * bladescalar(ej)
end

"""
    multivectorSum(ei::AbstractGeometricAlgebraType, ej::AbstractGeometricAlgebraType)::AbstractGeometricAlgebraType

Function that sums two multivectors and return its result.

# Arguments
- `ei::AbstractGeometricAlgebraType` : A multivector.
- `ej::AbstractGeometricAlgebraType` : A multivector.

# Return
The result Multivector.

"""
function multivectorSum(ei::AbstractGeometricAlgebraType, ej::AbstractGeometricAlgebraType)::AbstractGeometricAlgebraType
    result = ei.val + ej.val
    return Multivectors(result.nzind, result.nzval)
end


"""
    multivectorSub(ei::AbstractGeometricAlgebraType, ej::AbstractGeometricAlgebraType)::AbstractGeometricAlgebraType

Function that subtracts two multivectors and return its result.

# Arguments
- `ei::AbstractGeometricAlgebraType` : A multivector.
- `ej::AbstractGeometricAlgebraType` : A multivector.

# Return
The result Multivector.

"""
function multivectorSub(ei::AbstractGeometricAlgebraType, ej::AbstractGeometricAlgebraType)::AbstractGeometricAlgebraType
    result = ei.val - ej.val
    return Multivectors(result.nzind, result.nzval)
end

"""
    multivectorGP(ei::AbstractGeometricAlgebraType, ej::AbstractGeometricAlgebraType)::AbstractGeometricAlgebraType

Function that computes the geometric product of two multivectors and return its result.

# Arguments
- `ei::AbstractGeometricAlgebraType` : A multivector.
- `ej::AbstractGeometricAlgebraType` : A multivector.

# Return
The result Multivector.

"""
function multivectorGP(ei:: AbstractGeometricAlgebraType, ej::AbstractGeometricAlgebraType)::AbstractGeometricAlgebraType
    result = Multivectors([1],[0])
    for i in eachindex(ei.val.nzind)
        for j in eachindex(ej.val.nzind)
            result += bladeGeometricProduct(Multivectors([ei.val.nzind[i]],[ei.val.nzval[i]]), Multivectors([ej.val.nzind[j]],[ej.val.nzval[j]]))
        end
    end
    return result
end

"""
    multivectorIP(ei::AbstractGeometricAlgebraType, ej::AbstractGeometricAlgebraType)::AbstractGeometricAlgebraType

Function that computes the inner product of two multivectors and return its result.

# Arguments
- `ei::AbstractGeometricAlgebraType` : A multivector.
- `ej::AbstractGeometricAlgebraType` : A multivector.

# Return
The result Multivector.

"""
function multivectorIP(ei:: AbstractGeometricAlgebraType, ej::AbstractGeometricAlgebraType)::AbstractGeometricAlgebraType
    result = Multivectors([1],[0])
    for i in eachindex(ei.val.nzind)
        for j in eachindex(ej.val.nzind)
            result += bladeInnerProduct(Multivectors([ei.val.nzind[i]],[ei.val.nzval[i]]), Multivectors([ej.val.nzind[j]],[ej.val.nzval[j]]))
        end
    end
    return result
end

"""
    multivectorOP(ei::AbstractGeometricAlgebraType, ej::AbstractGeometricAlgebraType)::AbstractGeometricAlgebraType

Function that computes the outer product of two multivectors and return its result.

# Arguments
- `ei::AbstractGeometricAlgebraType` : A multivector.
- `ej::AbstractGeometricAlgebraType` : A multivector.

# Return
The result Multivector.

"""
function multivectorOP(ei:: AbstractGeometricAlgebraType, ej::AbstractGeometricAlgebraType)::AbstractGeometricAlgebraType
    result = Multivectors([1],[0])
    for i in eachindex(ei.val.nzind)
        for j in eachindex(ej.val.nzind)
            result += bladeOuterProduct(Multivectors([ei.val.nzind[i]],[ei.val.nzval[i]]), Multivectors([ej.val.nzind[j]],[ej.val.nzval[j]]))
        end
    end
    return result
end

"""
    multivectorByScalar(ei:: AbstractGeometricAlgebraType, k::Number)::AbstractGeometricAlgebraType

Function that computes the product by scalar of a multivectors and a scalar number and return its result.

# Arguments
- `ei::AbstractGeometricAlgebraType` : A multivector.
- `k::Number` : A number.

# Return
The result Multivector.

"""
function multivectorByScalar(ei:: AbstractGeometricAlgebraType, k::Number)::AbstractGeometricAlgebraType
    ei.val.nzval .*= k
    return ei
end