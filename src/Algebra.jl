include("CanonicalBasis.jl")

abstract type Algebra end

"""
    AlgebraFull(p, q, r, symbols, basis, basis_bit_order, metric, max)

A structure to define an algebra to be worked with its respective dimensions and canonical vectors.

# Fields
- `p::Int` : Represents the ammount of positive dimensions
- `q::Int` : Represents the ammount of negative dimensions
- `r::Int` : Represents the ammount of zero dimensions
- `symbols::Vector{String}` : Array of primary symbols for the Algebra
- `basis::Vector{String}` : Array of all symbols for the Algebra, normal order
- `basis_bit_order::Vector{String}` : Array of symbols for the Algebra, bit order
- `metric::Vector{Int8}` : Another way of representing the algebra signature
- `max::Integer` : max number of Algebra elements, the same as 2^(p+q+r)

"""
struct AlgebraFull{T <: Integer} <: Algebra
    p::Int
    q::Int
    r::Int
    symbols::Vector{String}
    basis::Vector{String}
    basis_bit_order::Vector{String}
    metric::Vector{Int8}
    max::T
end

"""
    AlgebraMin(p, q, r, symbols, metric, max)

A structure to define an algebra to be worked with its respective dimensions and canonical vectors.

# Fields
- `p::Int` : Represents the ammount of positive dimensions
- `q::Int` : Represents the ammount of negative dimensions
- `r::Int` : Represents the ammount of zero dimensions
- `symbols::Vector{String}` : Array of primary symbols for the Algebra
- `metric::Vector{Int8}` : Another way of representing the algebra signature
- `max::Int` : max number of Algebra elements, the same as 2^(p+q+r)

"""
struct AlgebraMin{T <: Integer} <: Algebra
    p::Int
    q::Int
    r::Int
    symbols::Vector{String}
    metric::Vector{Int8}
    max::T
end

"""
    describe(al::Algebra)

Describe function for showing the Algebra function.

# Arguments
- `al::Algebra` : The algebra for printing

"""
function describe(al::AlgebraFull)
    println("Algebra:")
    println("- p: $(al.p)")
    println("- q: $(al.q)")
    println("- r: $(al.r)")
    println("- symbols: $(al.symbols)")
    println("- basis: $(al.basis)")
    println("- basis_bit_order: $(al.basis_bit_order)")
    println("- metric: $(al.metric)")
    print("- max: $(al.max)")
end

function describe(al::AlgebraMin)
    println("Algebra:")
    println("- p: $(al.p)")
    println("- q: $(al.q)")
    println("- r: $(al.r)")
    println("- symbols: $(al.symbols)")
    println("- metric: $(al.metric)")
    print("- max: $(al.max)")
end

function Base.show(io::IO, al::AlgebraFull)
    describe(al)
end

function Base.show(io::IO, al::AlgebraMin)
    println("Algebra:")
    println("- p: $(al.p)")
    println("- q: $(al.q)")
    print("- r: $(al.r)")
end

"""
    create_algebra(p, [q], [r], [symbols])

Constructor Function of an algebraic object with signature p, q, r. If not defined, 
the symbols for the algebra are automatically calculated as canonical.

# Arguments
- `p::Int` : Represents the ammount of positive dimensions
- `q::Int` : Represents the ammount of negative dimensions
- `r::Int` : Represents the ammount of zero dimensions
- `symbols::Vector{String}` : Array of primary symbols for the Algebra

# Return
Returns the created Algebra object.

"""
function create_algebra(p::Int, q::Int = 0, r::Int = 0, symbols = nothing)::Algebra

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
        basis = sort(basis_bit_order, by = v -> length(v)) # May be removed for performance!!
    else
        if(length(symbols) != p+q+r)
            throw(DomainError(symbols, "The parameter 'symbols' has an incorrect length (should be equal to p+q+r)"))
        end
        basis_bit_order = canon_basis_bit_order(symbols)
        basis = canon_basis(symbols) # May be removed for performance!!
    end

    # Zero first for convention
    metric::Array{Int8} = vcat(fill(0, r), fill(1, p), fill(-1, q))
    max::Int = 1<<(p+q+r) # Fine to be Int

    global gb_current_algebra = AlgebraFull(p, q, r, symbols, basis, basis_bit_order, metric, max)
    return gb_current_algebra
end

function create_algebra_min(p, q = 0, r = 0)::Algebra

    if(p < 0)
        throw(DomainError(p, "The parameter 'p' must be greater than or equal to 0"))
    elseif(q < 0)
        throw(DomainError(q, "The parameter 'q' must be greater than or equal to 0"))
    elseif(r < 0)
        throw(DomainError(r, "The parameter 'r' must be greater than or equal to 0"))
    end

    metric::Array{Int8} = vcat(fill(0, r), fill(1, p), fill(-1, q))

    if(p+q+r <= 50000)
        symbols = canon_symbols(p, q, r)
    else
        symbols = []
    end

    if(p+q+r > 62) # Limit number for Int64, switching to BigInt!!
        max = BigInt(1)<<(p+q+r)
    else
        max = Int64(1)<<(p+q+r)
    end

    global gb_current_algebra = AlgebraMin{typeof(max)}(p, q, r, symbols, metric, max)
    return gb_current_algebra
end

# Global Variable for exporting the current Algebra
create_algebra(0)