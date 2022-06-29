angles = range(-2 * π, 2 * π, step=π/18)
angle_pairs = zip(range(-π, π, step=π/18), range(2π, 0, step=-π/18))

angle_good_complexes = map(θ -> Complex(cos(θ), sin(θ)), angles)
other_good_complexes = (Complex(1), im, Complex(-1), -im)
good_complexes = (angle_good_complexes..., other_good_complexes...)
bad_complexes = (Complex(0), Complex(2), Complex(-1.5), 2 * im, 1 + im)
complexes = (good_complexes..., bad_complexes...)

good_rotmats = map(θ -> [cos(θ) -sin(θ); sin(θ) cos(θ)], angles)
bad_rotmats = (
    [1 1; 0 1],
    zeros(2, 2),
    zeros(3, 3),
    ones(2, 2),
    -ones(2, 2))

@testset "constructor" begin
    for c = complexes
        @test SO2(c).c == c
    end
end

@testset "so2_from_complex" begin
    for c = complexes
        @test so2_from_complex(c, checks=false).c == c
    end

    for c = good_complexes
        @test so2_from_complex(c).c == c
    end

    for c = bad_complexes
        @test_throws DomainError so2_from_complex(c)
    end
end

@testset "so2_from_angle" begin
    for θ = angles
        @test so2_from_angle(θ).c == Complex(cos(θ), sin(θ))
    end
end

@testset "so2_from_rotmat" begin
    for r = (good_rotmats..., bad_rotmats...)
        @test so2_from_rotmat(r, checks=false).c == Complex(r[1, 1], r[2, 1])
    end

    for r = good_rotmats
        @test so2_from_rotmat(r).c == Complex(r[1, 1], r[2, 1])
    end

    for r = bad_rotmats
        @test_throws DomainError so2_from_rotmat(r)
    end
end

@testset "one" begin
    for θ = angles
        @test one(so2_from_angle(θ)).c == one(Complex)
    end

    @test one(SO2).c == one(Complex)
end

@testset "complex" begin
    for c = good_complexes
        @test complex(SO2(c)) == c
    end
end

@testset "angle" begin
    θ2c(θ) = Complex(cos(θ), sin(θ))
    for θ = angles
        @test θ2c(angle(so2_from_angle(θ))) ≈ θ2c(θ)
    end
end

@testset "rotmat" begin
    for r = good_rotmats
        @test rotmat(so2_from_rotmat(r)) == r
    end
end

@testset "multiplication" begin
    for θ = angles
        @test complex(one(SO2) * so2_from_angle(θ)) == Complex(cos(θ), sin(θ))
        @test complex(so2_from_angle(θ) * one(SO2)) == Complex(cos(θ), sin(θ))
    end

    for ap = angle_pairs
        θ1, θ2 = ap
        θs = θ1 + θ2
        l1, l2 = so2_from_angle(θ1), so2_from_angle(θ2)
        @test complex(l1 * l2) ≈ Complex(cos(θs), sin(θs))
    end
end

@testset "inverse" begin
    for θ = angles
        l = so2_from_angle(θ)
        @test complex(inv(l)) == Complex(cos(θ), -sin(θ))
        @test complex(inv(l) * l) ≈ complex(one(SO2))
        @test complex(l * inv(l)) ≈ complex(one(SO2))
    end
end

@testset "exponential map" begin
    for θ = angles
        @test complex(exp(SO2, fill(θ))) == Complex(cos(θ), sin(θ))
    end
end

@testset "logarithm" begin
    function test_log(l, θ)
        θl = log(l)[]
        @test θl <= π
        @test θl >= -π
        @test isapprox(sin(θl), sin(θ), atol=1e-15)
        @test isapprox(cos(θl), cos(θ), atol=1e-15)
    end

    for θ = angles 
        l = so2_from_angle(θ)
        test_log(l, θ)
        test_log(exp(SO2, fill(θ)), θ)
        @test complex(exp(SO2, log(l))) ≈ Complex(cos(θ), sin(θ))
    end
end
