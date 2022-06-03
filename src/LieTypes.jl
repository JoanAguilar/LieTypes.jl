module LieTypes

# Types
export LieGroup
export LieScalar

# Contsructors
export from_number

# Selectors
export number

# Operators include '*', 'inv', 'exp', and 'log', from Base.

include("abstract.jl")
include("scalar.jl")

end
