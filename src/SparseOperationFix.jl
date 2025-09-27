# Some fix in the Julia SparseArrays functions for BigInt
# I paste one of the functions and fixed the issue in this file

using SparseArrays

function sparsevec_operate( f::Function, 
                            x::SparseArrays.SparseVector{Tx, Tz}, 
                            y::SparseArrays.SparseVector{Ty, Tz}) where {Tx, Ty, Tz}

    n = length(x)
    xnzind = SparseArrays.nonzeroinds(x)
    xnzval = SparseArrays.nonzeros(x)
    ynzind = SparseArrays.nonzeroinds(y)
    ynzval = SparseArrays.nonzeros(y)
    R = Base.Broadcast.combine_eltypes(f, (x, y))
    mx = length(xnzind)
    my = length(ynzind)
    rind = Vector{Integer}(undef, mx+my)
    rval = Vector{R}(undef, mx+my)

    Base.require_one_based_indexing(xnzind, ynzind, xnzval, ynzval, rind, rval)

    ir = 0; ix = 1; iy = 1
    @inbounds while ix <= mx && iy <= my
        jx = xnzind[ix]
        jy = ynzind[iy]
        if jx == jy
            v = f(xnzval[ix], ynzval[iy])
            if SparseArrays._isnotzero(v)
                ir += 1; rind[ir] = jx; rval[ir] = v
            end
            ix += 1; iy += 1
        elseif jx < jy
            v = f(xnzval[ix], zero(Tx))
            ir += 1; rind[ir] = jx; rval[ir] = v
            ix += 1
        else
            v = f(zero(Tx), ynzval[iy])
            ir += 1; rind[ir] = jy; rval[ir] = v
            iy += 1
        end
    end
    @inbounds while ix <= mx
        v = f(xnzval[ix], zero(Ty))
        ir += 1; rind[ir] = xnzind[ix]; rval[ir] = v
        ix += 1
    end
    @inbounds while iy <= my
        v = f(zero(Ty), ynzval[iy])
        ir += 1; rind[ir] = ynzind[iy]; rval[ir] = v
        iy += 1
    end

    resize!(rind, ir)
    resize!(rval, ir)

    return sparsevec(rind, rval, n)

end