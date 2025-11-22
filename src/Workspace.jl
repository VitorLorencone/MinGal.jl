include("OperatorOverloading.jl")

"""
    change_algebra(al::Algebra)

Change the algebra REPL symbols on use. Not recomended because of time complexity.
Prefer using dot notation when working with multiple algebras.

# Arguments
- `al::Algebra`

"""
function change_algebra(al::Algebra)
    #TODO nova documentação e TESTES
    if typeof(al.max) == BigInt
        if al.p+ al.q + al.r >= 50000
            create_special_symbols(al)
        else
            create_symbols_min(al)
        end
    else
        create_symbols(al)
    end
    
    global gb_current_algebra = al
end

"""
    create_special_symbols()

Create and add to REPL the special symbols for this Algebra, such as
the "id" for scalar GAType grade 0 and "eI" for the pseudoscalar of the space. 

"""
function create_special_symbols(al::Algebra)

    mv1 = Multivector([0], [1], al)
    mv2 = Multivector([al.max-1], [1], al)

    al.blades[:id] = mv1
    al.blades[:eI] = mv2

    eval(:(id = $mv1))
    eval(:(export id))
    eval(:(eI = $mv2))
    eval(:(export eI))

end

"""
    create_symbols(stringSymbols)

Create and add to REPL all the (custom or basis) symbols for this Algebra.

# Arguments
- `string_symbols::Vector{String}` : An array with all the custom or basis symbols.

"""
function create_symbols(al::Algebra)

    symbol_array = []
    string_symbols = al.basis_bit_order

    for i in eachindex(string_symbols)
        if i == 1
            continue
        end
        push!(symbol_array, Symbol(string_symbols[i]))
    end

    for k in eachindex(symbol_array)
        symbol = symbol_array[k]
        mv = Multivector([k], [1], al)
        al.blades[symbol] = mv
        eval(:($symbol = $mv))
        eval(:(export $symbol))
    end

    create_special_symbols(al)
end

"""
    create_symbols_min(stringSymbols)

Create and add to REPL the (custom or basis) symbols for this Algebra. But
only the canonical basis blades symbols no combination between then.

# Arguments
- `string_symbols::Vector{String}` : An array with all the custom or basis symbols.

"""
function create_symbols_min(al::Algebra)

    symbol_array = []
    string_symbols = al.symbols

    for i in eachindex(string_symbols)
        push!(symbol_array, Symbol(string_symbols[i]))
    end

    for k in eachindex(symbol_array)
        symbol = symbol_array[k]

        if typeof(al.max) == BigInt
            mv = Multivector([BigInt(1)<<(k-1)], [1], al)
            al.blades[symbol] = mv
            eval(:($symbol = $mv))
        else
            mv = Multivector([1<<(k-1)], [1], al)
            al.blades[symbol] = mv
            eval(:($symbol = $mv))
        end
        
        eval(:(export $symbol))
    end

    create_special_symbols(al)
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
            al = create_algebra(p, q, r, symbols)
            create_symbols(al)

        elseif type == "min"
            al = create_algebra_min(p, q, r)
            create_symbols_min(al)

        elseif type == "special"
            al = create_algebra_min(p, q, r)
            create_special_symbols(al)

        else
            throw(Error(type, "Type does not exist."))
        end
        
    else
        if p+q+r >= 15
            al = create_algebra_min(p, q, r)
            
            if(p+q+r <= 50000)
                create_symbols_min(al)
            else
                create_special_symbols(al)
            end

        else
            al = create_algebra(p, q, r, symbols)
            create_symbols(al)
        end
    end

    return al
end