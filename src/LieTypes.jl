module LieTypes

# Types
export LieGroup
export LieScalar, LieVector

# Contsructors
export from_number, from_array

# Selectors
export number, array

# Operators include '*', 'inv', 'exp', and 'log', from Base.

include("abstract.jl")
include("scalar.jl")
include("vector.jl")

end
