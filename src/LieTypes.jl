module LieTypes

import LinearAlgebra as LA

# Support types and methods
export DualComplex, Quaternion
export scalar, abcd, sijk, vector, dual, iscomplex

# Types
export LieGroup
export LieScalar, LieVector, SO2, SE2, SO3

# Contsructors
export
    lie_scalar_from_number,
    lie_vector_from_array,
    so2_from_complex,
    so2_from_angle,
    so2_from_rotmat,
    se2_from_dual_complex,
    se2_from_matrix,
    se2_from_so2_disp,
    so3_from_quaternion,
    so3_from_rotmat,
    so3_from_axis_angle,
    so3_from_rotvec
# The use of standard constructors is discouraged.

# Selectors
export
    number,
    array,
    rotmat,
    dual_complex,
    matrix,
    so2,
    disp,
    quaternion,
    axis,
    rotvec
# As well as 'complex' and 'angle', from Base.

# Operators include '*', 'inv', 'exp', and 'log', from Base.

include("numbers/dualcomplex.jl")
include("numbers/quaternion.jl")
include("abstract.jl")
include("scalar.jl")
include("vector.jl")
include("so2.jl")
include("se2.jl")
include("so3.jl")

end
