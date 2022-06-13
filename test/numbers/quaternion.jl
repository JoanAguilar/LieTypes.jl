function test_values(q, s, i, j, k; comp=(==))
    @test comp(q.s, s)
    @test comp(q.i, i)
    @test comp(q.j, j)
    @test comp(q.k, k)
end

@testset "constructor" begin
    s, i, j, k = 1.0, 2.0, 3.0, 4.0

    test_values(Quaternion{Float32}(s, i, j, k), s, i, j, k)
    @test typeof(Quaternion{Float32}(s, i, j, k)) == Quaternion{Float32}

    test_values(Quaternion(s, i, j, k), s, i, j, k)
    test_values(Quaternion(s), s, 0.0, 0.0, 0.0)
    test_values(Quaternion(s, i=i), s, i, 0.0, 0.0)
    test_values(Quaternion(s, j=j), s, 0.0, j, 0.0)
    test_values(Quaternion(s, k=k), s, 0.0, 0.0, k)

    test_values(Quaternion(Complex(s, i)), s, i, 0.0, 0.0)
    test_values(Quaternion(Complex{Float32}(s, i)), s, i, 0.0, 0.0)
    @test typeof(Quaternion(Complex{Float32}(s, i))) == Quaternion{Float32} 
end

@testset "zero" begin
    test_values(zero(Quaternion{Float32}), 0.0, 0.0, 0.0, 0.0)
    @test typeof(zero(Quaternion{Float32})) == Quaternion{Float32}
    test_values(zero(Quaternion), 0, 0, 0, 0)
    test_values(zero(Quaternion(2)), 0, 0, 0, 0)
end

@testset "one" begin
    test_values(one(Quaternion{Float32}), 1.0, 0.0, 0.0, 0.0)
    @test typeof(one(Quaternion{Float32})) == Quaternion{Float32}
    test_values(one(Quaternion), 1, 0, 0, 0)
    test_values(one(Quaternion(2)), 1, 0, 0, 0)
end

@testset "scalar" begin
    @test scalar(zero(Quaternion)) == 0
    @test scalar(one(Quaternion)) == 1
    @test scalar(Quaternion(2, 3, 4, 5)) == 2
end

@testset "real" begin
    @test real(zero(Quaternion)) == 0
    @test real(one(Quaternion)) == 1
    @test real(Quaternion(2, 3, 4, 5)) == 2
end

@testset "complex" begin
    @test complex(zero(Quaternion)) == Complex(0, 0)
    @test complex(one(Quaternion)) == Complex(1, 0)
    @test complex(Quaternion(2, 3, 4, 5)) == Complex(2, 3)
end

@testset "sijk" begin
    @test sijk(zero(Quaternion)) == (0, 0, 0, 0)
    @test sijk(one(Quaternion)) == (1, 0, 0, 0)
    @test sijk(Quaternion(2, 3, 4, 5)) == (2, 3, 4, 5)
end

@testset "imag" begin
    @test imag(zero(Quaternion)) == 0
    @test imag(one(Quaternion)) == 0
    @test imag(Quaternion(2, 3, 4, 5)) == 3
end

@testset "vector" begin
    @test vector(zero(Quaternion)) == [0, 0, 0]
    @test vector(one(Quaternion)) == [0, 0, 0]
    @test vector(Quaternion(2, 3, 4, 5)) == [3, 4, 5]
end

@testset "isreal" begin
    @test isreal(zero(Quaternion))
    @test isreal(one(Quaternion))
    @test !isreal(Quaternion(2, 3, 0, 0))
    @test !isreal(Quaternion(2, 0, 4, 0))
    @test !isreal(Quaternion(2, 0, 0, 5))
    @test !isreal(Quaternion(2, 3, 4, 5))
end

@testset "iscomplex" begin
    @test iscomplex(zero(Quaternion))
    @test iscomplex(one(Quaternion))
    @test iscomplex(Quaternion(2, 3, 0, 0))
    @test !iscomplex(Quaternion(2, 0, 4, 0))
    @test !iscomplex(Quaternion(2, 0, 0, 5))
    @test !iscomplex(Quaternion(2, 3, 4, 5))
end

@testset "abs2" begin
    @test abs2(zero(Quaternion)) == 0
    @test abs2(one(Quaternion)) == 1
    @test abs2(Quaternion(2, 3, 4, 5)) == 54
end

@testset "abs" begin
    @test abs(zero(Quaternion)) == 0
    @test abs(one(Quaternion)) == 1
    @test abs(Quaternion(2, 3, 4, 5)) == √54
end

@testset "conj" begin
    test_values(conj(zero(Quaternion)), 0, 0, 0, 0)
    test_values(conj(one(Quaternion)), 1, 0, 0, 0)
    test_values(conj(Quaternion(2, 3, 4, 5)), 2, -3, -4, -5)
end

@testset "inv" begin
    test_values(inv(zero(Quaternion)), NaN, NaN, NaN, NaN, comp=isequal)
    test_values(inv(one(Quaternion)), 1, 0, 0, 0)
    test_values(inv(Quaternion(2, 3, 4, 5)), 2/54, -3/54, -4/54, -5/54)
end

@testset "addition" begin
    test_values(
        Quaternion(1, 2, 3, 4) + Quaternion(5, 6, 7, 8),
        6,
        8,
        10,
        12)
    test_values(1 + Quaternion(2, 3, 4, 5), 3, 3, 4, 5)
    test_values(Quaternion(1, 2, 3, 4) + 5, 6, 2, 3, 4)
    test_values(Complex(1, 2) + Quaternion(3, 4, 5, 6), 4, 6, 5, 6)
    test_values(Quaternion(1, 2, 3, 4) + Complex(5, 6), 6, 8, 3, 4)
    test_values(Quaternion(1, 2, 3, 4) + zero(Quaternion), 1, 2, 3, 4)
    test_values(zero(Quaternion) + Quaternion(1, 2, 3, 4), 1, 2, 3, 4)
end

@testset "subtraction" begin
    test_values(-Quaternion(1, 2, 3, 4), -1, -2, -3, -4)
    test_values(
        Quaternion(1, 2, 3, 4) - Quaternion(5, 6, 7, 8),
        -4,
        -4,
        -4,
        -4)
    test_values(1 - Quaternion(2, 3, 4, 5), -1, -3, -4, -5)
    test_values(Quaternion(1, 2, 3, 4) - 5, -4, 2, 3, 4)
    test_values(Complex(1, 2) - Quaternion(3, 4, 5, 6), -2, -2, -5, -6)
    test_values(Quaternion(1, 2, 3, 4) - Complex(5, 6), -4, -4, 3, 4)
    test_values(Quaternion(1, 2, 3, 4) - zero(Quaternion), 1, 2, 3, 4)
    test_values(zero(Quaternion) - Quaternion(1, 2, 3, 4), -1, -2, -3, -4)
end

@testset "multiplication" begin
    test_values(
        Quaternion(1, 2, 3, 4) * Quaternion(5, 6, 7, 8),
        -60,
        12,
        30,
        24)
    test_values(1 * Quaternion(2, 3, 4, 5), 2, 3, 4, 5)
    test_values(Quaternion(1, 2, 3, 4) * 5, 5, 10, 15, 20)
    test_values(Complex(1, 2) * Quaternion(3, 4, 5, 6), -5, 10, -7, 16)
    test_values(Quaternion(1, 2, 3, 4) * Complex(5, 6), -7, 16, 39, 2)
    test_values(Quaternion(1, 2, 3, 4) * zero(Quaternion), 0, 0, 0, 0)
    test_values(zero(Quaternion) * Quaternion(1, 2, 3, 4), 0, 0, 0, 0)
    test_values(Quaternion(1, 2, 3, 4) * one(Quaternion), 1, 2, 3, 4)
    test_values(one(Quaternion) * Quaternion(1, 2, 3, 4), 1, 2, 3, 4)
end

@testset "division" begin
    @test Quaternion(1, 2, 3, 4) / Quaternion(5, 6, 7, 8) ≈
        Quaternion(1, 2, 3, 4) * inv(Quaternion(5, 6, 7, 8))
    @test 1 / Quaternion(2, 3, 4, 5) ≈ inv(Quaternion(2, 3, 4, 5))
    test_values(Quaternion(1, 2, 3, 4) / 5, 1/5, 2/5, 3/5, 4/5)
    @test Complex(1, 2) / Quaternion(3, 4, 5, 6) ≈
        Complex(1, 2) * inv(Quaternion(3, 4, 5, 6))
    @test Quaternion(1, 2, 3, 4) / Complex(5, 6) ≈
        Quaternion(1, 2, 3, 4) * inv(Quaternion(Complex(5, 6)))
    test_values(
        Quaternion(1, 2, 3, 4) / zero(Quaternion),
        NaN,
        NaN,
        NaN,
        NaN,
        comp=isequal)
    test_values(zero(Quaternion) / Quaternion(1, 2, 3, 4), 0, 0, 0, 0)
    test_values(Quaternion(1, 2, 3, 4) / one(Quaternion), 1, 2, 3, 4)
    @test one(Quaternion) / Quaternion(1, 2, 3, 4) ≈
        inv(Quaternion(1, 2, 3, 4))

    @test Quaternion(1, 2, 3, 4) \ Quaternion(5, 6, 7, 8) ≈
        Quaternion(5, 6, 7, 8) / Quaternion(1, 2, 3, 4)
    @test 1 \ Quaternion(2, 3, 4, 5) ≈ Quaternion(2, 3, 4, 5)
    @test Quaternion(1, 2, 3, 4) \ 5 ≈ 5 / Quaternion(1, 2, 3, 4)
    @test Complex(1, 2) \ Quaternion(3, 4, 5, 6) ≈
        Quaternion(3, 4, 5, 6) / Complex(1, 2)
    @test Quaternion(1, 2, 3, 4) \ Complex(5, 6) ≈
        Complex(5, 6) / Quaternion(1, 2, 3, 4)
    @test Quaternion(1, 2, 3, 4) \ zero(Quaternion) == zero(Quaternion)
    test_values(
        zero(Quaternion) \ Quaternion(1, 2, 3, 4),
        NaN,
        NaN,
        NaN,
        NaN,
        comp=isequal)
    @test Quaternion(1, 2, 3, 4) \ one(Quaternion) ≈
        inv(Quaternion(1, 2, 3, 4))
    @test one(Quaternion) \ Quaternion(1, 2, 3, 4) == Quaternion(1, 2, 3, 4)
end

@testset "==" begin
    @test Quaternion(1, 2, 3, 4) == Quaternion(1, 2, 3, 4)
    @test Quaternion(1, 2, 3, 4) != Quaternion(5, 2, 3, 4)
    @test Quaternion(1, 2, 3, 4) != Quaternion(1, 5, 3, 4)
    @test Quaternion(1, 2, 3, 4) != Quaternion(1, 2, 5, 4)
    @test Quaternion(1, 2, 3, 4) != Quaternion(1, 2, 3, 5)
end

@testset "norm" begin
    @test LA.norm(Quaternion(1, 2, 3, 4)) == abs(Quaternion(1, 2, 3, 4))
    @test LA.norm(Quaternion(1, 2, 3, 4), 1) == abs(Quaternion(1, 2, 3, 4))
    @test LA.norm(Quaternion(1, 2, 3, 4), 3) == abs(Quaternion(1, 2, 3, 4))
end

@testset "normalize" begin
    @test LA.norm(LA.normalize(Quaternion(1, 2, 3, 4))) ≈ 1
    @test LA.norm(LA.normalize(Quaternion(5, 6, 7, 8))) ≈ 1
end
