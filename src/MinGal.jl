# TODO
# -> Implement the multiplicative approach
# -> More GA Functions (Dual)
# -> Fix Docs and Tests
# -> Linear Systems Examples

module MinGal

include("Workspace.jl")

export Algebra, grade, grade_projection, reverse
export geometric_product, inner_product, outer_product
export bitmap, scalar, multivector_sum, multivector_sub
export Multivector, Blade, product_by_scalar, describe

end