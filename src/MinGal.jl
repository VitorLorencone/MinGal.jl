# TODO
# -> Implement the multiplicative approach
# -> Redo Operation Table
# -> More GA Functions
# -> Fix creating huge dimensions
# -> Fix Docs and Tests
# -> Per-grade compression and new data structure
# -> Linear Systems Examples

module MinGal

include("Workspace.jl")

export Algebra, grade, grade_projection
export geometric_product, inner_product, outer_product
export bitmap, scalar, multivector_sum, multivector_sub
export Multivector, Blade, product_by_scalar, describe

end