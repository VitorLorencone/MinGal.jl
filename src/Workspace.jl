include("OperatorOverloading.jl")

# Create "id" value for any algebra
const id = Multivector([0], [1])

"""
    create_symbols(stringSymbols)

Create and add to REPL all the (custom or basis) symbols for this Algebra.

# Arguments
- `string_symbols::Array` : An array with all the custom or basis symbols.

"""
function create_symbols(string_symbols::Array)

    symbol_array = []

    for i in eachindex(string_symbols)
        if i == 1
            continue
        end
        push!(symbol_array, Symbol(string_symbols[i]))
    end

    for k in eachindex(symbol_array)
        symbol = symbol_array[k]
        eval(:($symbol = Multivector([$k], [1])))
        eval(:(export $symbol))
    end

end

"""
    Algebra(p, q, r, symbols)::Algebra

Main function for creating your Algebra and adding its basis blades to REPL.
Constructor Function of an algebraic object with parameters p, q, R^{p, q}, and its multivector space.
If not defined, the last two parameters are automatically calculated as canonical.

# Arguments
- `p::Int` : The first parameter of the definition
- `q::Int` : The second parameter of the definition
- `r::Int` : The third parameter of the definition
- `symbols::Array{String}` : An Array with string vectors to work with

# Return
Returns the created Algebra object.

"""
function Algebra(p = 0, q = 0, r = 0, symbols = nothing)::Algebra

    Al = create_algebra(p, q, r, symbols)
    create_symbols(Al.basis_bit_order)

    return Al

end