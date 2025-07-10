# TODO
# -> Implement the multiplicative approach
# -> Redo Operation Table
# -> More GA Functions
# -> Fix creating huge dimensions
# -> Dinamically create symbols

module MinGal

include("Workspace.jl")

export Algebra, id, grade, grade_projection
export geometric_product, inner_product, outer_product
export bitmap, scalar, multivector_sum, multivector_sub
export Multivector, Blade, product_by_scalar, describe

end