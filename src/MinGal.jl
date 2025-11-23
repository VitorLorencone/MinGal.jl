# TODO
# -> Implement the multiplicative approach
# -> Optimize some Geometric Functions
# -> More examples

module MinGal

include("Workspace.jl")

export Algebra, grade, grade_projection, revert
export geometric_product, inner_product, outer_product
export bitmap, scalar, multivector_sum, multivector_sub
export Multivector, Blade, product_by_scalar, describe
export grade_selection, left_contraction, right_contraction
export GAType, GAVector, GAArray, grade_involution, invert, dual
export undual, regressive_product, exp_ga, is_scalar, is_blade
export canonical_basis, chain, norm, reverse_norm, euclidean_norm
export conjugate, clifford_conjugation, conjugate_norm, clifford_norm
export grade_minus, grade_plus, grade_null, scalar_product
export euclidean_scalar_product, change_algebra

end