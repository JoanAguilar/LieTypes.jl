arrays = (
    fill(0),
    fill(-0.5),
    [1],
    [-1.5, 2.0],
    [-2 3],
    [-2.5 3.5],
    [-3 4; 5 -6],
    [-3.5 4.5; 5.5 -6.5],
    ones(5, 5))

pairs = (
    (fill(-0.5), fill(1.75)),
    ([-1.75], [0.5]),
    ([1.0, 1.5], [-2.5, -3.75]),
    ([1.0 -1.5; -0.5 -2.0], [2.5 0.5; 3.0 3.5]))

@testset "constructor" begin
    for a = arrays
        @test LieVector(a).a == a 
    end
end

@testset "lie_vector_from_array" begin
    for a = arrays
        @test lie_vector_from_array(a).a == a
    end
end

@testset "one" begin
    for a = arrays
        @test one(LieVector(a)).a == zero(a)
    end
end

@testset "array" begin
    for a = arrays
        @test array(LieVector(a)) == a
    end
end

@testset "multiplication" begin
    for a = arrays
        @test (one(LieVector(a)) * LieVector(a)).a == a
        @test (LieVector(a) * one(LieVector(a))).a == a
    end

    for p = pairs
        a1, a2 = p
        l1, l2 = LieVector(a1), LieVector(a2)
        @test (l1 * l2).a == a1 + a2
    end
end

@testset "inverse" begin
    for a = arrays
        l = LieVector(a)
        @test inv(l).a == -a
        @test (inv(l) * l).a == zero(a)
        @test (l * inv(l)).a == zero(a)
    end
end

@testset "exponential map" begin
    for a = arrays
        @test exp(LieVector, a).a == a
    end
end

@testset "logarithm" begin
    for a = arrays
        l = LieVector(a)
        @test log(l) == a
        @test log(exp(LieVector, a)) == a
        @test exp(LieVector, log(l)).a == a
    end
end
