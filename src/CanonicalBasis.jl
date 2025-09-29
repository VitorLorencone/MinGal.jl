import Combinatorics
using .Combinatorics

"""
    canon_symbols(p, [q], [r])::Vector{String}

Function that writes the canonical vector space symbols, given the parameters p, q and r for definition

# Arguments
- `p::Int` : Represents the ammount of positive dimensions
- `q::Int` : Represents the ammount of negative dimensions
- `r::Int` : Represents the ammount of zero dimensions

# Return
Return an array of strings with all the necessary elements for this space.

"""
function canon_symbols(p::Int, q::Int = 0, r::Int = 0)::Vector{String}
    
    if(p < 0)
        throw(DomainError(p,"The parameter 'p' must be greater than or equal to 0"))
    elseif(q < 0)
        throw(DomainError(q,"The parameter 'q' must be greater than or equal to 0"))
    elseif(r < 0)
        throw(DomainError(r,"The parameter 'r' must be greater than or equal to 0"))
    end

    if(r == 0)
        starting_index = 1
    else
        starting_index = 0
    end

    basis::Vector{String} = []

    for i in starting_index:(p+q+r-1+starting_index)
        str_name::String = "e" * string(i)
        push!(basis, str_name)
    end

    return basis
end

"""
    canon_basis(symbols)::Vector{String}

Function that lists all the combinations of canonical vectors in a given Algebra.

# Arguments
- `symbols::Vector{String}` : An array of strings to be combined.

# Return
Returns a list with all combinations of the elements, forming the basis of the multivector space.

"""
function canon_basis(symbols::Vector{String})::Vector{String}
    
    basis::Array{String} = ["1"]
    for k in 1:length(symbols)
        for val in combinations(symbols, k)
            push!(basis, join(val))
        end
    end

    return basis
end

"""
    canon_basis_bit_order(symbols)::Vector{String}

Function that lists all the combinations of canonical vectors in a given Algebra in bit order.

# Arguments
- `symbols::Vector{String}` : An array of strings to be combined.

# Return
Returns a list with all combinations of the elements, forming the basis of the multivector space, in bit order.

"""
function canon_basis_bit_order(symbols)::Vector{String}
    n = length(symbols)
    basis = ["1"]
    for i in 1:(1<<n - 1)
        val = String[]
        for j in 0:(n-1)
            if (i >> j) & 1 == 1
                push!(val, symbols[j+1])
            end
        end
        push!(basis, join(val))
    end
    return basis
end

# TODO this is useless and inneficient by now. The best approach is to
# create the map while creating one of the canon basis
"""
    binary_index_map(base, basis)::Dict{Int, Int}

Function that maps the index of normal and bit order.

# Arguments
- `base::Vector{String}` : An array of strings to be combined.
- `basis::Vector{String}` : An array of the final combination of strings.

# Return
Returns a dict with all the mapped values.

"""
function binary_index_map(base::Vector{String}, basis::Vector{String})
    n = length(base)
    map1 = Dict{Int, Int}()
    map2 = Dict{Int, Int}()
    map1[0] = 1
    map2[1] = 0
    for i in 1:(1<<n - 1)
        bits = digits(i, base=2, pad=n)
        str = join([base[j] for j in 1:n if bits[j] == 1])
        k = findfirst(==(str), basis)
        map1[i] = k
        map2[k] = i
    end

    global gb_bin_to_nat = map1
    global gb_nat_to_bin = map2
end