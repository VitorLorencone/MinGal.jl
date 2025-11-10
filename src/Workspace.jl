include("OperatorOverloading.jl")

function create_special_symbols()
    eval(:(id = Multivector([0], [1])))
    eval(:(export id))
    eval(:(eI = Multivector([gb_current_algebra.max-1], [1])))
    eval(:(export eI))
end

"""
    create_symbols(stringSymbols)

Create and add to REPL all the (custom or basis) symbols for this Algebra.

# Arguments
- `string_symbols::Vector{String}` : An array with all the custom or basis symbols.

"""
function create_symbols(string_symbols::Vector{String})

    symbol_array = []

    for i in eachindex(string_symbols)
        if i == 1
            continue
        end
        push!(symbol_array, Symbol(string_symbols[i]))
    end

    for k in eachindex(symbol_array)
        symbol = symbol_array[k]
        mv = Multivector([k], [1])
        eval(:($symbol = $mv))
        eval(:(export $symbol))
    end

    create_special_symbols()
end

function create_symbols_min(string_symbols::Vector{String})

    symbol_array = []

    for i in eachindex(string_symbols)
        push!(symbol_array, Symbol(string_symbols[i]))
    end

    for k in eachindex(symbol_array)
        symbol = symbol_array[k]

        if typeof(gb_current_algebra.max) == BigInt
            mv = Multivector([BigInt(1)<<(k-1)], [1])
            eval(:($symbol = $mv))
        else
            mv = Multivector([1<<(k-1)], [1])
            eval(:($symbol = $mv))
        end
        
        eval(:(export $symbol))
    end

    create_special_symbols()
end

"""
    Algebra(p, q, r, symbols)::Algebra

Main function for creating your Algebra and adding its basis blades to REPL.
Constructor Function of an algebraic object with parameters p, q, r, R^{p, q, r}, and its multivector space.
If not defined, the last parameter is automatically calculated as canonical.

# Arguments
- `p::Int` : The first parameter of the definition
- `q::Int` : The second parameter of the definition
- `r::Int` : The third parameter of the definition
- `symbols::Vector{String}` : An Array with string vectors to work with

# Return
Returns the created Algebra object.

"""
function Algebra(p = 0, q = 0, r = 0, symbols = nothing, type = nothing)::Algebra

    if type !== nothing

        if type == "full"
            Al = create_algebra(p, q, r, symbols)
            create_symbols(Al.basis_bit_order)

        elseif type == "min"
            Al = create_algebra_min(p, q, r)
            create_symbols_min(Al.symbols)

        elseif type == "special"
            Al = create_algebra_min(p, q, r)
            create_special_symbols()

        else
            throw(Error(type, "Type does not exist."))
        end
        
    else
        if p+q+r >= 15
            Al = create_algebra_min(p, q, r)
            
            if(p+q+r <= 50000)
                create_symbols_min(Al.symbols)
            else
                create_special_symbols()
            end

        else
            Al = create_algebra(p, q, r, symbols)
            create_symbols(Al.basis_bit_order)
        end
    end

    return Al
end