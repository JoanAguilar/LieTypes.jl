module LieTypes

import LinearAlgebra as LA

# Types
export LieGroup
export LieScalar, LieVector, SO2

# Contsructors
export from_number, from_array, from_complex, from_angle, from_rotmat

# Selectors
export number, array, rotmat
# As well as 'complex' and 'angle', from Base.

# Operators include '*', 'inv', 'exp', and 'log', from Base.

include("abstract.jl")
include("scalar.jl")
include("vector.jl")
include("so2.jl")

end
