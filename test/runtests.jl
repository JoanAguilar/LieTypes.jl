using LieTypes
using Test

import LinearAlgebra as LA

@testset "LieTypes.jl" begin
    @testset "DualComplex" begin include("numbers/dualcomplex.jl") end
    @testset "LieScalar" begin include("scalar.jl") end
    @testset "LieVector" begin include("vector.jl") end
    @testset "SO2" begin include("so2.jl") end
    @testset "SE2" begin include("se2.jl") end
end
