function test_values(n, r, d; comp=(==))
    @test comp(n.r, r)
    @test comp(n.d, d)
end

@testset "constructor" begin
    r = 1.0
    d = 2.0

    test_values(Dual{Float32}(r, d), r, d)
    @test typeof(Dual{Float32}(r, d)) == Dual{Float32}

    test_values(Dual(r, d), r, d)
    test_values(Dual(r), r, 0.0)
    test_values(Dual(r, d=d), r, d)
end

@testset "zero" begin
    for T = (Float32, Complex{Float64}, Quaternion{Float64})
        test_values(zero(Dual{T}), zero(T), zero(T))
        @test typeof(zero(Dual{T})) == Dual{T}
    end
    test_values(zero(Dual), 0, 0)
    test_values(zero(Dual(2)), 0, 0)
end

@testset "one" begin
    for T = (Float32, Complex{Float64}, Quaternion{Float64})
        test_values(one(Dual{T}), one(T), zero(T))
        @test typeof(one(Dual{T})) == Dual{T}
    end
    test_values(one(Dual), 1, 0)
    test_values(one(Dual(2)), 1, 0)
end

@testset "dual unit" begin
    test_values(du, 0, 1)
    test_values(ϵ, 0, 1)
end

@testset "real" begin
    @test real(zero(Dual)) == 0
    @test real(one(Dual)) == 1
    @test real(Dual(2, 3)) == 2
    @test real(Dual(Complex(1, 2), Complex(3, 4))) == Complex(1, 2)
    @test real(Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8))) ==
        Quaternion(1, 2, 3, 4)
end

@testset "true real" begin
    @test truereal(zero(Dual)) == 0
    @test truereal(one(Dual)) == 1
    @test truereal(Dual(2, 3)) == 2
    @test truereal(Dual(Complex(1, 2), Complex(3, 4))) == 1
    @test truereal(Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8))) == 1
end

@testset "redu" begin
    @test redu(zero(Dual)) == (0, 0)
    @test redu(one(Dual)) == (1, 0)
    @test redu(Dual(2, 3)) == (2, 3)
    @test redu(Dual(Complex(1, 2), Complex(3, 4))) ==
        (Complex(1, 2), Complex(3, 4))
    @test redu(Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8))) ==
        (Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8))
end

@testset "dual" begin
    @test dual(zero(Dual)) == 0
    @test dual(one(Dual)) == 0
    @test dual(Dual(2, 3)) == 3
    @test dual(Dual(Complex(1, 2), Complex(3, 4))) == Complex(3, 4)
    @test dual(Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8))) ==
        Quaternion(5, 6, 7, 8)
end

@testset "isreal" begin
    @test isreal(zero(Dual))
    @test isreal(one(Dual))
    @test isreal(Dual(Complex(1)))
    @test isreal(Dual(Quaternion(1)))
    @test !isreal(Dual(2, 3))
    @test !isreal(Dual(0, 3))
    @test !isreal(Dual(Complex(0, 1)))
    @test !isreal(Dual(Quaternion(0, 1, 0, 0)))
    @test !isreal(Dual(Quaternion(0, 0, 1, 0)))
    @test !isreal(Dual(Quaternion(0, 0, 0, 1)))
end

@testset "abs2" begin
    @test abs2(zero(Dual)) == 0
    @test abs2(one(Dual)) == 1
    @test abs2(Dual(2, 3)) == Dual(4, 12)
    @test abs2(Dual(Complex(1, 2), Complex(3, 4))) == Dual(5, 22)
    @test (==)(
        abs2(Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8))),
        Dual(30, 140))
end

@testset "conj" begin
    test_values(conj(zero(Dual)), 0, 0)
    test_values(conj(one(Dual)), 1, 0)
    test_values(conj(Dual(2, 3)), 2, 3)
    test_values(
        conj(Dual(Complex(1, 2), Complex(3, 4))),
        Complex(1, -2),
        Complex(3, -4))
    test_values(
        conj(Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8))),
        Quaternion(1, -2, -3, -4),
        Quaternion(5, -6, -7, -8))
end

@testset "dual conj" begin
    test_values(dualconj(zero(Dual)), 0, 0)
    test_values(dualconj(one(Dual)), 1, 0)
    test_values(dualconj(Dual(2, 3)), 2, -3)
    test_values(
        dualconj(Dual(Complex(1, 2), Complex(3, 4))),
        Complex(1, 2),
        Complex(-3, -4))
    test_values(
        dualconj(Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8))),
        Quaternion(1, 2, 3, 4),
        Quaternion(-5, -6, -7, -8))
end


@testset "inv" begin
    test_values(inv(zero(Dual)), Inf, NaN, comp=isequal)
    test_values(inv(one(Dual)), 1, 0)
    test_values(inv(Dual(2, 3)), 2/4, -3/4)
    test_values(
        inv(Dual(Complex(1, 2), Complex(3, 4))),
        Complex(0.2, -0.4),
        Complex(-0.28, 0.96),
        comp=isapprox)
    test_values(
        inv(Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8))),
        Quaternion(1 / 30, -2 / 30, -3 / 30, -4 / 30),
        Quaternion(1 / 90, 1 / 9, 7 / 30, 16/ 45),
        comp=isapprox)
end

@testset "addition" begin
    test_values(Dual(1, 2) + Dual(3, 4), 4, 6)
    test_values(
        Dual(Complex(1, 2), Complex(3, 4)) +
            Dual(Complex(5, 6), Complex(7, 8)),
        Complex(6, 8),
        Complex(10, 12))
    test_values(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) +
            Dual(Quaternion(9, 10, 11, 12), Quaternion(13, 14, 15, 16)),
        Quaternion(10, 12, 14, 16),
        Quaternion(18, 20, 22, 24))
    test_values(1 + Dual(2, 3), 3, 3)
    test_values(Dual(1, 2) + 3, 4, 2)
    test_values(
        1 + Dual(Complex(2, 3), Complex(4, 5)),
        Complex(3, 3),
        Complex(4, 5))
    test_values(
        Dual(Complex(1, 2), Complex(3, 4)) + 5,
        Complex(6, 2),
        Complex(3, 4))
    test_values(
        1 + Dual(Quaternion(2, 3, 4, 5), Quaternion(6, 7, 8, 9)),
        Quaternion(3, 3, 4, 5),
        Quaternion(6, 7, 8, 9))
    test_values(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) + 9,
        Quaternion(10, 2, 3, 4),
        Quaternion(5, 6, 7, 8))
    test_values(ϵ + Dual(2, 3), 2, 4)
    test_values(Dual(1, 2) + 3ϵ, 1, 5)
    test_values(
        ϵ + Dual(Complex(2, 3), Complex(4, 5)),
        Complex(2, 3),
        Complex(5, 5))
    test_values(
        Dual(Complex(1, 2), Complex(3, 4)) + 5ϵ,
        Complex(1, 2),
        Complex(8, 4))
    test_values(
        ϵ + Dual(Quaternion(2, 3, 4, 5), Quaternion(6, 7, 8, 9)),
        Quaternion(2, 3, 4, 5),
        Quaternion(7, 7, 8, 9))
    test_values(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) + 9ϵ,
        Quaternion(1, 2, 3, 4),
        Quaternion(14, 6, 7, 8))
    test_values(Dual(1, 2) + zero(Dual), 1, 2)
    test_values(zero(Dual) + Dual(1, 2), 1, 2)
    test_values(
        Dual(Complex(1, 2), Complex(3, 4)) + zero(Dual{Complex}),
        Complex(1, 2),
        Complex(3, 4))
    test_values(
        zero(Dual{Complex}) + Dual(Complex(1, 2), Complex(3, 4)),
        Complex(1, 2),
        Complex(3, 4))
    test_values(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) +
            zero(Dual{Quaternion}),
        Quaternion(1, 2, 3, 4),
        Quaternion(5, 6, 7, 8))
    test_values(
        zero(Dual{Quaternion}) +
            Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)),
        Quaternion(1, 2, 3, 4),
        Quaternion(5, 6, 7, 8))
end

@testset "subtraction" begin
    test_values(-Dual(1, 2), -1, -2)
    test_values(
        Dual(Complex(1, 2), Complex(3, 4)) -
            Dual(Complex(5, 6), Complex(7, 8)),
        Complex(-4, -4),
        Complex(-4, -4))
    test_values(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) -
            Dual(Quaternion(9, 10, 11, 12), Quaternion(13, 14, 15, 16)),
        Quaternion(-8, -8, -8, -8),
        Quaternion(-8, -8, -8, -8))
    test_values(
        -Dual(Complex(1, 2), Complex(3, 4)),
        Complex(-1, -2),
        Complex(-3, -4))
    test_values(
        -Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)),
        Quaternion(-1, -2, -3, -4),
        Quaternion(-5, -6, -7, -8))

    test_values(Dual(1, 2) - Dual(3, 4), -2, -2)
    test_values(1 - Dual(2, 3), -1, -3)
    test_values(Dual(1, 2) - 3, -2, 2)
    test_values(
        1 - Dual(Complex(2, 3), Complex(4, 5)),
        Complex(-1, -3),
        Complex(-4, -5))
    test_values(
        Dual(Complex(1, 2), Complex(3, 4)) - 5,
        Complex(-4, 2),
        Complex(3, 4))
    test_values(
        1 - Dual(Quaternion(2, 3, 4, 5), Quaternion(6, 7, 8, 9)),
        Quaternion(-1, -3, -4, -5),
        Quaternion(-6, -7, -8, -9))
    test_values(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) - 9,
        Quaternion(-8, 2, 3, 4),
        Quaternion(5, 6, 7, 8))
    test_values(ϵ - Dual(2, 3), -2, -2)
    test_values(Dual(1, 2) - 3ϵ, 1, -1)
    test_values(
        ϵ - Dual(Complex(2, 3), Complex(4, 5)),
        Complex(-2, -3),
        Complex(-3, -5))
    test_values(
        Dual(Complex(1, 2), Complex(3, 4)) - 5ϵ,
        Complex(1, 2),
        Complex(-2, 4))
    test_values(
        ϵ - Dual(Quaternion(2, 3, 4, 5), Quaternion(6, 7, 8, 9)),
        Quaternion(-2, -3, -4, -5),
        Quaternion(-5, -7, -8, -9))
    test_values(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) - 9ϵ,
        Quaternion(1, 2, 3, 4),
        Quaternion(-4, 6, 7, 8))
    test_values(Dual(1, 2) - zero(Dual), 1, 2)
    test_values(zero(Dual) - Dual(1, 2), -1, -2)
    test_values(
        Dual(Complex(1, 2), Complex(3, 4)) - zero(Dual{Complex}),
        Complex(1, 2),
        Complex(3, 4))
    test_values(
        zero(Dual{Complex}) - Dual(Complex(1, 2), Complex(3, 4)),
        Complex(-1, -2),
        Complex(-3, -4))
    test_values(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) -
            zero(Dual{Quaternion}),
        Quaternion(1, 2, 3, 4),
        Quaternion(5, 6, 7, 8))
    test_values(
        zero(Dual{Quaternion}) -
            Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)),
        Quaternion(-1, -2, -3, -4),
        Quaternion(-5, -6, -7, -8))
end

@testset "multiplication" begin
    test_values(Dual(1, 2) * Dual(3, 4), 3, 10)
    test_values(
        Dual(Complex(1, 2), Complex(3, 4)) *
            Dual(Complex(5, 6), Complex(7, 8)),
        Complex(-7, 16),
        Complex(-18, 60))
    test_values(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) *
            Dual(Quaternion(9, 10, 11, 12), Quaternion(13, 14, 15, 16)),
        Quaternion(-92, 20, 54, 40),
        Quaternion(-312, 128, 204, 184))
    test_values(1 * Dual(2, 3), 2, 3)
    test_values(Dual(1, 2) * 3, 3, 6)
    test_values(
        1 * Dual(Complex(2, 3), Complex(4, 5)),
        Complex(2, 3),
        Complex(4, 5))
    test_values(
        Dual(Complex(1, 2), Complex(3, 4)) * 5,
        Complex(5, 10),
        Complex(15, 20))
    test_values(
        1 * Dual(Quaternion(2, 3, 4, 5), Quaternion(6, 7, 8, 9)),
        Quaternion(2, 3, 4, 5),
        Quaternion(6, 7, 8, 9))
    test_values(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) * 9,
        Quaternion(9, 18, 27, 36),
        Quaternion(45, 54, 63, 72))
    test_values(ϵ * Dual(2, 3), 0, 2)
    test_values(Dual(1, 2) * 3ϵ, 0, 3)
    test_values(
        ϵ * Dual(Complex(2, 3), Complex(4, 5)),
        Complex(0, 0),
        Complex(2, 3))
    test_values(
        Dual(Complex(1, 2), Complex(3, 4)) * 5ϵ,
        Complex(0, 0),
        Complex(5, 10))
    test_values(
        ϵ * Dual(Quaternion(2, 3, 4, 5), Quaternion(6, 7, 8, 9)),
        Quaternion(0, 0, 0, 0),
        Quaternion(2, 3, 4, 5))
    test_values(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) * 9ϵ,
        Quaternion(0, 0, 0, 0),
        Quaternion(9, 18, 27, 36))
    test_values(Dual(1, 2) * one(Dual), 1, 2)
    test_values(one(Dual) * Dual(1, 2), 1, 2)
    test_values(
        Dual(Complex(1, 2), Complex(3, 4)) * one(Dual{Complex}),
        Complex(1, 2),
        Complex(3, 4))
    test_values(
        one(Dual{Complex}) * Dual(Complex(1, 2), Complex(3, 4)),
        Complex(1, 2),
        Complex(3, 4))
    test_values(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) *
            one(Dual{Quaternion}),
        Quaternion(1, 2, 3, 4),
        Quaternion(5, 6, 7, 8))
    test_values(
        one(Dual{Quaternion}) *
            Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)),
        Quaternion(1, 2, 3, 4),
        Quaternion(5, 6, 7, 8))
end

@testset "division" begin
    @test Dual(1, 2) / Dual(3, 4) ≈ Dual(1, 2) * inv(Dual(3, 4))
    @test isapprox(
        Dual(Complex(1, 2), Complex(3, 4)) /
            Dual(Complex(5, 6), Complex(7, 8)),
        Dual(Complex(1, 2), Complex(3, 4)) *
            inv(Dual(Complex(5, 6), Complex(7, 8))))
    @test isapprox(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) /
            Dual(Quaternion(9, 10, 11, 12), Quaternion(13, 14, 15, 16)),
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) *
            inv(Dual(Quaternion(9, 10, 11, 12), Quaternion(13, 14, 15, 16))))
    @test 1 / Dual(2, 3) ≈ inv(Dual(2, 3))
    @test isapprox(
        1 / Dual(Complex(2, 3), Complex(4, 5)),
        inv(Dual(Complex(2, 3), Complex(4, 5))))
    @test isapprox(
        1 / Dual(Quaternion(2, 3, 4, 5), Quaternion(6, 7, 8, 9)),
        inv(Dual(Quaternion(2, 3, 4, 5), Quaternion(6, 7, 8, 9))))
    test_values(Dual(1, 2) / 3, 1/3, 2/3)
    @test isapprox(
        Dual(Complex(1, 2), Complex(3, 4)) / 5,
        Dual(Complex(1 / 5, 2 / 5), Complex(3 / 5, 4 / 5)))
    @test isapprox(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) / 9,
        Dual(
            Quaternion(1 / 9, 2 / 9, 3 / 9, 4 /9),
            Quaternion(5 / 9, 6 / 9, 7 / 9, 8 / 9)))
    test_values(Dual(1, 2) / zero(Dual), Inf, NaN, comp=isequal)
    test_values(
        Dual(Complex(1, 2), Complex(3, 4)) / zero(Dual{Complex}),
        Complex(NaN, NaN),
        Complex(NaN, NaN),
        comp=isequal)
    test_values(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) /
            zero(Dual{Quaternion}),
        Quaternion(NaN, NaN, NaN, NaN),
        Quaternion(NaN, NaN, NaN, NaN),
        comp=isequal)
    test_values(zero(Dual) / Dual(1, 2), 0, 0)
    test_values(zero(Dual{Complex}) / Dual(Complex(1, 2), Complex(3, 4)), 0, 0)
    test_values(
        zero(Dual{Quaternion}) /
            Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)),
        0,
        0)
    test_values(Dual(1, 2) / one(Dual), 1, 2)
    test_values(
        Dual(Complex(1, 2), Complex(3, 4)) / one(Dual{Complex}),
        Complex(1, 2),
        Complex(3, 4))
    test_values(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) /
            one(Dual{Quaternion}),
        Quaternion(1, 2, 3, 4),
        Quaternion(5, 6, 7, 8))
    @test one(Dual) / Dual(1, 2) ≈ inv(Dual(1, 2))
    @test isapprox(
        one(Dual{Complex}) / Dual(Complex(1, 2), Complex(3, 4)),
        inv(Dual(Complex(1, 2), Complex(3, 4))))
    @test isapprox(
        one(Dual{Quaternion}) /
            Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)),
        inv(Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8))))

    @test Dual(1, 2) \ Dual(3, 4) ≈ inv(Dual(1, 2)) * Dual(3, 4)
    @test isapprox(
        Dual(Complex(1, 2), Complex(3, 4)) \
            Dual(Complex(5, 6), Complex(7, 8)),
        inv(Dual(Complex(1, 2), Complex(3, 4))) *
            Dual(Complex(5, 6), Complex(7, 8)))
    @test isapprox(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) \
            Dual(Quaternion(9, 10, 11, 12), Quaternion(13, 14, 15, 16)),
        inv(Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8))) *
            Dual(Quaternion(9, 10, 11, 12), Quaternion(13, 14, 15, 16)))
    test_values(1 \ Dual(2, 3), 2, 3)
    test_values(
        1 \ Dual(Complex(2, 3), Complex(4, 5)),
        Complex(2, 3),
        Complex(4, 5))
    test_values(
        1 \ Dual(Quaternion(2, 3, 4, 5), Quaternion(6, 7, 8, 9)),
        Quaternion(2, 3, 4, 5),
        Quaternion(6, 7, 8, 9))
    @test Dual(1, 2) \ 3 ≈ inv(Dual(1, 2)) * 3
    @test isapprox(
        Dual(Complex(1, 2), Complex(3, 4)) \ 5,
        inv(Dual(Complex(1, 2), Complex(3, 4))) * 5)
    @test isapprox(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) \ 9,
        inv(Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8))) * 9)
    test_values(Dual(1, 2) \ zero(Dual), 0, 0)
    test_values(Dual(Complex(1, 2), Complex(3, 4)) \ zero(Dual{Complex}), 0, 0)
    test_values(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) \
            zero(Dual{Quaternion}),
        0,
        0)
    test_values(zero(Dual) \ Dual(1, 2), Inf, NaN, comp=isequal)
    test_values(
        zero(Dual{Complex}) \ Dual(Complex(1, 2), Complex(3, 4)),
        Complex(NaN, NaN),
        Complex(NaN, NaN),
        comp=isequal)
    test_values(
        zero(Dual{Quaternion}) \
            Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)),
        Quaternion(NaN, NaN, NaN, NaN),
        Quaternion(NaN, NaN, NaN, NaN),
        comp=isequal)
    @test Dual(1, 2) \ one(Dual) ≈ inv(Dual(1, 2))
    @test isapprox(
        Dual(Complex(1, 2), Complex(3, 4)) \ one(Dual{Complex}),
        inv(Dual(Complex(1, 2), Complex(3, 4))))
    @test isapprox(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)) \
            one(Dual{Quaternion}),
        inv(Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8))))
    test_values(one(Dual) \ Dual(1, 2), 1, 2)
    test_values(
        one(Dual{Complex}) \ Dual(Complex(1, 2), Complex(3, 4)),
        Complex(1, 2),
        Complex(3, 4))
    test_values(
        one(Dual{Quaternion}) \
            Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)),
        Quaternion(1, 2, 3, 4),
        Quaternion(5, 6, 7, 8))
end

@testset "==" begin
    @test Dual(1, 2) == Dual(1, 2)
    @test (==)(
        Dual(Complex(1, 2), Complex(3, 4)),
        Dual(Complex(1, 2), Complex(3, 4)))
    @test (==)(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)),
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)))
    @test Dual(1, 2) != Dual(3, 2)
    @test Dual(1, 2) != Dual(1, 3)
    @test (!=)(
        Dual(Complex(1, 2), Complex(3, 4)),
        Dual(Complex(5, 2), Complex(3, 4)))
    @test (!=)(
        Dual(Complex(1, 2), Complex(3, 4)),
        Dual(Complex(1, 5), Complex(3, 4)))
    @test (!=)(
        Dual(Complex(1, 2), Complex(3, 4)),
        Dual(Complex(1, 2), Complex(5, 4)))
    @test (!=)(
        Dual(Complex(1, 2), Complex(3, 4)),
        Dual(Complex(1, 2), Complex(3, 5)))
    @test (!=)(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)),
        Dual(Quaternion(9, 2, 3, 4), Quaternion(5, 6, 7, 8)))
    @test (!=)(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)),
        Dual(Quaternion(1, 2, 3, 4), Quaternion(9, 6, 7, 8)))
end

@testset "isapprox" begin
    @test isapprox(Dual(1, 2), Dual(1, 2))
    @test isapprox(
        Dual(Complex(1, 2), Complex(3, 4)),
        Dual(Complex(1, 2), Complex(3, 4)))
    @test isapprox(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)),
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)))
    @test !isapprox(Dual(1, 2), Dual(1.1, 2.1))
    @test !isapprox(
        Dual(Complex(1, 2), Complex(3, 4)),
        Dual(Complex(1.1, 2.1), Complex(3.1, 4.1)))
    @test !isapprox(
        Dual(Quaternion(1, 2, 3, 4), Quaternion(5, 6, 7, 8)),
        Dual(Quaternion(1.1, 2.1, 3.1, 4.1), Quaternion(5.1, 6.1, 7.1, 8.1)))
end
