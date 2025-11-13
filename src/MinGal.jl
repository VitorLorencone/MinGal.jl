# TODO
# -> Implement the multiplicative approach
# -> Optimize some Geometric Functions
# -> Update Docs and Tests
# -> More examples

module MinGal

include("Workspace.jl")

export Algebra, grade, grade_projection, revert
export geometric_product, inner_product, outer_product
export bitmap, scalar, multivector_sum, multivector_sub
export Multivector, Blade, product_by_scalar, describe
export grade_selection, left_contraction, right_contraction
export GAType, GAVector, GAArray, involution, invert, dual
export undual, regressive_product, exp_ga, isScalar, isBlade
export canonical_basis, chain, norm, reverse_norm, euclidian_norm

end