function test_values(dc, a, b, c, d; comp=(==))
    @test comp(dc.a, a)
    @test comp(dc.b, b)
    @test comp(dc.c, c)
    @test comp(dc.d, d)
end

@testset "constructor" begin
    a, b, c, d = 1.0, 2.0, 3.0, 4.0

    test_values(DualComplex{Float32}(a, b, c, d), a, b, c, d)
    @test typeof(DualComplex{Float32}(a, b, c, d)) == DualComplex{Float32}

    test_values(DualComplex(a, b, c, d), a, b, c, d)
    test_values(DualComplex(a), a, 0.0, 0.0, 0.0)
    test_values(DualComplex(a, b=b), a, b, 0.0, 0.0)
    test_values(DualComplex(a, c=c), a, 0.0, c, 0.0)
    test_values(DualComplex(a, d=d), a, 0.0, 0.0, d)

    test_values(DualComplex(Complex(a, b)), a, b, 0.0, 0.0)
    test_values(DualComplex(Complex{Float32}(a, b)), a, b, 0.0, 0.0)
    @test typeof(DualComplex(Complex{Float32}(a, b))) == DualComplex{Float32} 
end

@testset "zero" begin
    test_values(zero(DualComplex{Float32}), 0.0, 0.0, 0.0, 0.0)
    @test typeof(zero(DualComplex{Float32})) == DualComplex{Float32}
    test_values(zero(DualComplex), 0, 0, 0, 0)
    test_values(zero(DualComplex(2)), 0, 0, 0, 0)
end

@testset "one" begin
    test_values(one(DualComplex{Float32}), 1.0, 0.0, 0.0, 0.0)
    @test typeof(one(DualComplex{Float32})) == DualComplex{Float32}
    test_values(one(DualComplex), 1, 0, 0, 0)
    test_values(one(DualComplex(2)), 1, 0, 0, 0)
end

@testset "scalar" begin
    @test scalar(zero(DualComplex)) == 0
    @test scalar(one(DualComplex)) == 1
    @test scalar(DualComplex(2, 3, 4, 5)) == 2
end

@testset "real" begin
    @test real(zero(DualComplex)) == 0
    @test real(one(DualComplex)) == 1
    @test real(DualComplex(2, 3, 4, 5)) == 2
end

@testset "complex" begin
    @test complex(zero(DualComplex)) == Complex(0, 0)
    @test complex(one(DualComplex)) == Complex(1, 0)
    @test complex(DualComplex(2, 3, 4, 5)) == Complex(2, 3)
end

@testset "abcd" begin
    @test abcd(zero(DualComplex)) == (0, 0, 0, 0)
    @test abcd(one(DualComplex)) == (1, 0, 0, 0)
    @test abcd(DualComplex(2, 3, 4, 5)) == (2, 3, 4, 5)
end

@testset "imag" begin
    @test imag(zero(DualComplex)) == 0
    @test imag(one(DualComplex)) == 0
    @test imag(DualComplex(2, 3, 4, 5)) == 3
end

@testset "vector" begin
    @test vector(zero(DualComplex)) == [0, 0, 0]
    @test vector(one(DualComplex)) == [0, 0, 0]
    @test vector(DualComplex(2, 3, 4, 5)) == [3, 4, 5]
end

@testset "dual" begin
    @test dual(zero(DualComplex)) == (0, 0)
    @test dual(one(DualComplex)) == (0, 0)
    @test dual(DualComplex(2, 3, 4, 5)) == (4, 5)
end

@testset "abs2" begin
    @test abs2(zero(DualComplex)) == 0
    @test abs2(one(DualComplex)) == 1
    @test abs2(DualComplex(2, 3, 4, 5)) == 13
end

@testset "abs" begin
    @test abs(zero(DualComplex)) == 0
    @test abs(one(DualComplex)) == 1
    @test abs(DualComplex(2, 3, 4, 5)) == √13
end

@testset "conj" begin
    test_values(conj(zero(DualComplex)), 0, 0, 0, 0)
    test_values(conj(one(DualComplex)), 1, 0, 0, 0)
    test_values(conj(DualComplex(2, 3, 4, 5)), 2, -3, -4, -5)
end

@testset "inv" begin
    test_values(inv(zero(DualComplex)), NaN, NaN, NaN, NaN, comp=isequal)
    test_values(inv(one(DualComplex)), 1, 0, 0, 0)
    test_values(inv(DualComplex(2, 3, 4, 5)), 2/13, -3/13, -4/13, -5/13)
end

@testset "isreal" begin
    @test isreal(zero(DualComplex))
    @test isreal(one(DualComplex))
    @test !isreal(DualComplex(2, 3, 0, 0))
    @test !isreal(DualComplex(2, 3, 4, 5))
end

@testset "iscomplex" begin
    @test iscomplex(zero(DualComplex))
    @test iscomplex(one(DualComplex))
    @test iscomplex(DualComplex(2, 3, 0, 0))
    @test !iscomplex(DualComplex(2, 3, 4, 5))
end

@testset "addition" begin
    test_values(
        DualComplex(1, 2, 3, 4) + DualComplex(5, 6, 7, 8),
        6,
        8,
        10,
        12)
    test_values(1 + DualComplex(2, 3, 4, 5), 3, 3, 4, 5)
    test_values(DualComplex(1, 2, 3, 4) + 5, 6, 2, 3, 4)
    test_values(Complex(1, 2) + DualComplex(3, 4, 5, 6), 4, 6, 5, 6)
    test_values(DualComplex(1, 2, 3, 4) + Complex(5, 6), 6, 8, 3, 4)
    test_values(DualComplex(1, 2, 3, 4) + zero(DualComplex), 1, 2, 3, 4)
    test_values(zero(DualComplex) + DualComplex(1, 2, 3, 4), 1, 2, 3, 4)
end

@testset "subtraction" begin
    test_values(-DualComplex(1, 2, 3, 4), -1, -2, -3, -4)
    test_values(
        DualComplex(1, 2, 3, 4) - DualComplex(5, 6, 7, 8),
        -4,
        -4,
        -4,
        -4)
    test_values(1 - DualComplex(2, 3, 4, 5), -1, -3, -4, -5)
    test_values(DualComplex(1, 2, 3, 4) - 5, -4, 2, 3, 4)
    test_values(Complex(1, 2) - DualComplex(3, 4, 5, 6), -2, -2, -5, -6)
    test_values(DualComplex(1, 2, 3, 4) - Complex(5, 6), -4, -4, 3, 4)
    test_values(DualComplex(1, 2, 3, 4) - zero(DualComplex), 1, 2, 3, 4)
    test_values(zero(DualComplex) - DualComplex(1, 2, 3, 4), -1, -2, -3, -4)
end

@testset "multiplication" begin
    test_values(
        DualComplex(1, 2, 3, 4) * DualComplex(5, 6, 7, 8),
        -7,
        16,
        30,
        24)
    test_values(1 * DualComplex(2, 3, 4, 5), 2, 3, 4, 5)
    test_values(DualComplex(1, 2, 3, 4) * 5, 5, 10, 15, 20)
    test_values(Complex(1, 2) * DualComplex(3, 4, 5, 6), -5, 10, -7, 16)
    test_values(DualComplex(1, 2, 3, 4) * Complex(5, 6), -7, 16, 39, 2)
    test_values(DualComplex(1, 2, 3, 4) * zero(DualComplex), 0, 0, 0, 0)
    test_values(zero(DualComplex) * DualComplex(1, 2, 3, 4), 0, 0, 0, 0)
    test_values(DualComplex(1, 2, 3, 4) * one(DualComplex), 1, 2, 3, 4)
    test_values(one(DualComplex) * DualComplex(1, 2, 3, 4), 1, 2, 3, 4)
end

@testset "division" begin
    @test DualComplex(1, 2, 3, 4) / DualComplex(5, 6, 7, 8) ≈
        DualComplex(1, 2, 3, 4) * inv(DualComplex(5, 6, 7, 8))
    @test 1 / DualComplex(2, 3, 4, 5) ≈ inv(DualComplex(2, 3, 4, 5))
    test_values(DualComplex(1, 2, 3, 4) / 5, 1/5, 2/5, 3/5, 4/5)
    @test Complex(1, 2) / DualComplex(3, 4, 5, 6) ≈
        Complex(1, 2) * inv(DualComplex(3, 4, 5, 6))
    @test DualComplex(1, 2, 3, 4) / Complex(5, 6) ≈
        DualComplex(1, 2, 3, 4) * inv(DualComplex(Complex(5, 6)))
    test_values(
        DualComplex(1, 2, 3, 4) / zero(DualComplex),
        NaN,
        NaN,
        NaN,
        NaN,
        comp=isequal)
    test_values(zero(DualComplex) / DualComplex(1, 2, 3, 4), 0, 0, 0, 0)
    test_values(DualComplex(1, 2, 3, 4) / one(DualComplex), 1, 2, 3, 4)
    @test one(DualComplex) / DualComplex(1, 2, 3, 4) ≈
        inv(DualComplex(1, 2, 3, 4))

    @test DualComplex(1, 2, 3, 4) \ DualComplex(5, 6, 7, 8) ≈
        DualComplex(5, 6, 7, 8) / DualComplex(1, 2, 3, 4)
    @test 1 \ DualComplex(2, 3, 4, 5) ≈ DualComplex(2, 3, 4, 5)
    @test DualComplex(1, 2, 3, 4) \ 5 ≈ 5 / DualComplex(1, 2, 3, 4)
    @test Complex(1, 2) \ DualComplex(3, 4, 5, 6) ≈
        DualComplex(3, 4, 5, 6) / Complex(1, 2)
    @test DualComplex(1, 2, 3, 4) \ Complex(5, 6) ≈
        Complex(5, 6) / DualComplex(1, 2, 3, 4)
    @test DualComplex(1, 2, 3, 4) \ zero(DualComplex) == zero(DualComplex)
    test_values(
        zero(DualComplex) \ DualComplex(1, 2, 3, 4),
        NaN,
        NaN,
        NaN,
        NaN,
        comp=isequal)
    @test DualComplex(1, 2, 3, 4) \ one(DualComplex) ≈
        inv(DualComplex(1, 2, 3, 4))
    @test one(DualComplex) \ DualComplex(1, 2, 3, 4) == DualComplex(1, 2, 3, 4)
end

@testset "==" begin
    @test DualComplex(1, 2, 3, 4) == DualComplex(1, 2, 3, 4)
    @test DualComplex(1, 2, 3, 4) != DualComplex(5, 2, 3, 4)
    @test DualComplex(1, 2, 3, 4) != DualComplex(1, 5, 3, 4)
    @test DualComplex(1, 2, 3, 4) != DualComplex(1, 2, 5, 4)
    @test DualComplex(1, 2, 3, 4) != DualComplex(1, 2, 3, 5)
end

@testset "norm" begin
    @test LA.norm(DualComplex(1, 2, 3, 4)) == abs(DualComplex(1, 2, 3, 4))
    @test LA.norm(DualComplex(1, 2, 3, 4), 1) == abs(DualComplex(1, 2, 3, 4))
    @test LA.norm(DualComplex(1, 2, 3, 4), 3) == abs(DualComplex(1, 2, 3, 4))
end

@testset "normalize" begin
    @test LA.norm(LA.normalize(DualComplex(1, 2, 3, 4))) ≈ 1
    @test LA.norm(LA.normalize(DualComplex(5, 6, 7, 8))) ≈ 1
end
