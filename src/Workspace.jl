include("OperatorOverloading.jl")

# Create "id" value for any algebra
const id = Multivectors([1], [1])

"""
    CreateSymbols(stringSymbols)

Create and add to REPL all the (custom or basis) symbols for this Algebra.

# Arguments
- `stringSymbols::Array` : An array with all the custom or basis symbols.

"""
function CreateSymbols(stringSymbols::Array)

    symbolArray = []

    for i in eachindex(stringSymbols)
        if i == 1
            continue
        end
        push!(symbolArray, Symbol(stringSymbols[i][1]))
    end

    for k in eachindex(symbolArray)
        symbol = symbolArray[k]
        eval(:(const $symbol = Multivectors([$k+1], [1])))
        eval(:(export $symbol))
    end

end

"""
    CreateTables()

Create and add to REPL all the custom operation tables for this Algebra.

"""
function CreateTables() 
    CreateGPTable()
    CreateIPTable()
    CreateOPTable()
end

"""
    Algebra(p, q, VectorBasis, Basis)::AlgebraStruct

Main function for creating your Algebra and adding its basis blades to REPL.
Constructor Function of an algebraic object with parameters p, q, R^{p, q}, and its multivector space.
If not defined, the last two parameters are automatically calculated as canonical.

# Arguments
- `p::Int` : The first parameter of the definition
- `q::Int` : The second parameter of the definition
- `VectorBasis::Array{String}` : An Array with vectors to work with
- `Basis::Array{Tuple{String,Int}}` : An Array with the multivector base and it's indexes

# Return
Returns the created Algebra object.

"""
function Algebra(p = 0, q = 0, VectorBasis = CanonVectorBasis(p, q), Basis = CanonBasis(VectorBasis))::AlgebraStruct

    Al = CreateAlgebra(p, q, VectorBasis, Basis)
    CreateSymbols(Basis)
    CreateTables()

    return Al

end