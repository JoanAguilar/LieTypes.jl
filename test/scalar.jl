reals = (0, 1, -1, -1.0, 0.0, 1.0)
complexes = (Complex(0), 1 * im, 1 + im, Complex(0.0),  1.0 * im, 1.0 + im)
numbers = (reals..., complexes...)

pairs = (
    (0, 0),
    (1, 0),
    (-1.0, 0.0),
    (-1.0, -1.0),
    (0, Complex(0)),
    (1 * im, 1),
    (1.0, 1.0 + im),
    (-1 * im, 1 * im))

@testset "constructor" begin
    for n = numbers
        @test LieScalar(n).s == n 
    end
end

@testset "from_number" begin
    for n = numbers
        @test from_number(LieScalar, n).s == n
    end
end

@testset "one" begin
    for n = numbers
        @test one(LieScalar(n)).s == 0
    end

    @test one(LieScalar).s == 0
    for T = (Int64, Float64, Complex)
        @test one(LieScalar{T}).s == 0
    end
end

@testset "number" begin
    for n = numbers
        @test number(LieScalar(n)) == n
    end
end

@testset "multiplication" begin
    for n = numbers
        @test (one(LieScalar(n)) * LieScalar(n)).s == n
        @test (LieScalar(n) * one(LieScalar(n))).s == n
        @test (one(LieScalar) * LieScalar(n)).s == n
        @test (LieScalar(n) * one(LieScalar)).s == n
    end

    for p = pairs
        n1, n2 = p
        l1, l2 = LieScalar(n1), LieScalar(n2)
        @test (l1 * l2).s == n1 + n2
    end
end

@testset "inverse" begin
    for n = numbers
        l = LieScalar(n)
        @test inv(l).s == -n
        @test (inv(l) * l).s == 0
        @test (l * inv(l)).s == 0
    end
end

@testset "exponential map" begin
    for n = numbers
        a = fill(n)
        @test exp(LieScalar, a).s == n
    end
end

@testset "logarithm" begin
    for n = numbers
        l = LieScalar(n)
        @test log(l)[] == n
        @test log(exp(LieScalar, fill(n)))[] == n
        @test exp(LieScalar, log(l)).s == n
    end
end
