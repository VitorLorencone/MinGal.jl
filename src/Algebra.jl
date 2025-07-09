include("CanonicalBasis.jl")

"""
    Algebra(p, q, r, symbols, basis, basis_bit_order, metric, max)

A structure to define an algebra to be worked with its respective dimensions and canonical vectors.

# Fields
- `p::Int` : Represents the ammount of positive dimensions
- `q::Int` : Represents the ammount of negative dimensions
- `z::Int` : Represents the ammount of zero dimensions
- `symbols::Array{String}` : Array of primary symbols for the Algebra
- `basis::Array{String}` : Array of all symbols for the Algebra, normal order
- `basis_bit_order::Array{String}` : Array of symbols for the Algebra, bit order
- `metric::Array{Int8}` : Another way of representing the algebra signature
- `max::Int` : max number of Algebra, the same as 2^(p+q+r)

"""
struct Algebra
    p::Int
    q::Int
    r::Int
    symbols::Array{String}
    basis::Array{String}
    basis_bit_order::Array{String}
    metric::Array{Int8}
    max::Int
end

"""
    describe(al::Algebra)

Describe function for showing the Algebra function.

# Arguments
- `al::Algebra` : The algebra for printing

"""
function describe(al::Algebra)
    println("Algebra:")
    println("- p: $(al.p)")
    println("- q: $(al.q)")
    println("- r: $(al.r)")
    println("- symbols: $(al.symbols)")
    println("- basis: $(al.basis)")
    println("- basis_bit_order: $(al.basis_bit_order)")
    println("- metric: $(al.metric)")
    println("- max: $(al.max)")
end

"""
    create_algebra(p, [q], [r], [symbols])

Constructor Function of an algebraic object with signature p, q, r. If not defined, 
the symbols for the algebra are automatically calculated as canonical.

# Arguments
- `p::Int` : Represents the ammount of positive dimensions
- `q::Int` : Represents the ammount of negative dimensions
- `r::Int` : Represents the ammount of zero dimensions
- `symbols::Array{String}` : Array of primary symbols for the Algebra

# Return
Returns the created Algebra object.

"""
function create_algebra(p, q = 0, r = 0, symbols = nothing)::Algebra

    if(p < 0)
        throw(DomainError(p, "The parameter 'p' must be greater than or equal to 0"))
    elseif(q < 0)
        throw(DomainError(q, "The parameter 'q' must be greater than or equal to 0"))
    elseif(r < 0)
        throw(DomainError(r, "The parameter 'r' must be greater than or equal to 0"))
    end

    if symbols === nothing
        symbols = canon_symbols(p, q, r)
        basis_bit_order = canon_basis_bit_order(symbols)
        basis = sort(basis_bit_order, by = v -> length(v)) # May be commented by now
    else
        if(length(symbols) != p+q+r)
            throw(DomainError(symbols, "The parameter 'symbols' has an incorrect length (should be equal to p+q+r)"))
        end
        basis_bit_order = canon_basis_bit_order(symbols)
        basis = canon_basis(symbols) # May be commented by now
    end

    metric::Array{Int} = vcat(fill(0, r), fill(1, p), fill(-1, q))
    max::Int = 2^(p+q+r)

    global gb_current_algebra = Algebra(p, q, r, symbols, basis, basis_bit_order, metric, max)
    return gb_current_algebra
end

# Global Variable for exporting the current Algebra
create_algebra(0)