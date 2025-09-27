# TODO
# -> Implement the multiplicative approach
# -> More GA Functions (Dual)
# -> Fix Docs and Tests
# -> Linear Systems Examples
# -> Optimize some Geometric Functions

module MinGal

include("Workspace.jl")

export Algebra, grade, grade_projection, reverse
export geometric_product, inner_product, outer_product
export bitmap, scalar, multivector_sum, multivector_sub
export Multivector, Blade, product_by_scalar, describe
export grade_selection, left_contraction, right_contraction
export GAType, GAVector, GAArray, involution, invert, dual
export undual, regressive_product, exp_ga, isScalar, isBlade

end