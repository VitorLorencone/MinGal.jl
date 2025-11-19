include("OperatorOverloading.jl")

"""
    change_algebra(al::Algebra)

Change with safety the algebra on use. Does not recreate the symbols.

# Arguments
- `al::Algebra`

"""
function change_algebra(al::Algebra)
    global gb_current_algebra = al
end

"""
    create_special_symbols()

Create and add to REPL the special symbols for this Algebra, such as
the "id" for scalar GAType grade 0 and "eI" for the pseudoscalar of the space. 

"""
function create_special_symbols()
    eval(:(id = Multivector([0], [1])))
    eval(:(export id))
    eval(:(eI = Multivector([gb_current_algebra.max-1], [1])))
    eval(:(export eI))
    gb_current_algebra.blades[:id] = Multivector([0], [1])
    gb_current_algebra.blades[:eI] = Multivector([gb_current_algebra.max-1], [1])
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
        gb_current_algebra.blades[symbol] = mv
        eval(:($symbol = $mv))
        eval(:(export $symbol))
    end

    create_special_symbols()
end

"""
    create_symbols_min(stringSymbols)

Create and add to REPL the (custom or basis) symbols for this Algebra. But
only the canonical basis blades symbols no combination between then.

# Arguments
- `string_symbols::Vector{String}` : An array with all the custom or basis symbols.

"""
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
            gb_current_algebra.blades[symbol] = mv
        else
            mv = Multivector([1<<(k-1)], [1])
            eval(:($symbol = $mv))
            gb_current_algebra.blades[symbol] = mv
        end
        
        eval(:(export $symbol))
    end

    create_special_symbols()
end

"""
    Algebra(p, [q], [r], [symbols], [type])::Algebra

Main function for creating your Algebra and adding its basis blades to REPL.
Constructor Function of an algebraic object with parameters p, q, r, R^{p, q, r}, and its multivector space.
If not defined, the last two parameters are automatically calculated as canonical.

# Arguments
- `p::Int` : Represents the ammount of positive dimensions
- `q::Int` : Represents the ammount of negative dimensions
- `r::Int` : Represents the ammount of null dimensions
- `symbols::Vector{String}` : Array of primary symbols for the Algebra
- `type::String` : String for manually selecting "full", "min" or "special" Algebras

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