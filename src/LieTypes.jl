module LieTypes

import LinearAlgebra as LA

# Support types and methods
export DualComplex, scalar, vector, dual, iscomplex

# Types
export LieGroup
export LieScalar, LieVector, SO2, SE2

# Contsructors
export
    from_number,
    from_array,
    from_complex,
    from_angle,
    from_rotmat,
    from_dual_complex,
    from_matrix,
    from_so2_disp
# The use of standard constructors is discouraged.

# Selectors
export number, array, rotmat, dual_complex, matrix, so2, disp
# As well as 'complex' and 'angle', from Base.

# Operators include '*', 'inv', 'exp', and 'log', from Base.

include("numbers/dualcomplex.jl")
include("abstract.jl")
include("scalar.jl")
include("vector.jl")
include("so2.jl")
include("se2.jl")

end
