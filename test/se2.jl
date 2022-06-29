function angle_disp_to_dual_complex(angle, disp)
    chθ, shθ = cos(angle/2), sin(angle/2)
    r = DualComplex(chθ, shθ, shθ, shθ)
    t = DualComplex(1, 0, -disp[2], disp[1])
    return t * r
end

function angle_disp_to_matrix(angle, disp)
    return [
        cos(angle) -sin(angle) disp[1];
        sin(angle) cos(angle) disp[2];
        0 0 1]
end

angles_disps = (
    (angle=0, disp=[0, 0]),
    (angle=π/6, disp=[0, 0]),
    (angle=0.75π, disp=[0, 0]),
    (angle=-π/6, disp=[0, 0]),
    (angle=-0.75π, disp=[0, 0]),
    (angle=0, disp=[1, 0]),
    (angle=0, disp=[0, 1]),
    (angle=0, disp=[1, 2]),
    (angle=π/6, disp=[-1, 2]),
    (angle=-π/3, disp=[3, -4]))

angles_disps_good_dual_complexes = map(
    ad -> angle_disp_to_dual_complex(ad.angle, ad.disp),
    angles_disps)
other_good_dual_complexes = (
    DualComplex(1),
    DualComplex(0, b=1),
    DualComplex(-1),
    DualComplex(0, b=-1),
    DualComplex(1, 0, 3, 4))
good_dual_complexes = (
    angles_disps_good_dual_complexes...,
    other_good_dual_complexes...)
bad_dual_complexes = (
    zero(DualComplex),
    DualComplex(2),
    DualComplex(-1.5),
    DualComplex(0, b=2),
    DualComplex(1, b=1),
    DualComplex(0, c=1, d=2))
dual_complexes = (good_dual_complexes..., bad_dual_complexes...)

good_matrices = map(
    ad -> angle_disp_to_matrix(ad.angle, ad.disp),
    angles_disps)
bad_matrices = (
    [1 1 2;
     0 1 3;
     0 0 1],
    [1 0 2;
     0 1 3;
     4 5 6],
    zeros(3, 3),
    ones(3, 3),
    zeros(4, 4))
matrices = (good_matrices..., bad_matrices...)

@testset "constructor" begin
    for dc = dual_complexes
        @test SE2(dc).dc == dc
    end
end

@testset "se2_from_dual_complex" begin
    for dc = dual_complexes
        @test se2_from_dual_complex(dc, checks=false).dc == dc
    end

    for dc = good_dual_complexes
        @test se2_from_dual_complex(dc).dc == dc
    end

    for dc = bad_dual_complexes
        @test_throws DomainError se2_from_dual_complex(dc)
    end
end

@testset "se2_from_matrix" begin
    for mat in matrices
        @test_nowarn se2_from_matrix(mat, checks=false)
    end

    for ad in angles_disps
        θ, d = ad
        dc = angle_disp_to_dual_complex(θ, d)
        mat = angle_disp_to_matrix(θ, d)
        @test abs(se2_from_matrix(mat).dc - dc) ≈ 0 atol=1e-16
    end

    for mat in bad_matrices
        @test_throws DomainError se2_from_matrix(mat)
    end
end

@testset "se2_from_so2_disp" begin
    for mat in matrices
        rotmat = mat[1:2, 1:2]
        r = so2_from_rotmat(rotmat, checks=false)
        disp = mat[1:2, 3]
        @test_nowarn se2_from_so2_disp(r, disp, checks=false)
    end

    for ad in angles_disps
        θ, disp = ad
        dc = angle_disp_to_dual_complex(θ, disp)
        mat = angle_disp_to_matrix(θ, disp)
        rotmat = mat[1:2, 1:2]
        r = so2_from_rotmat(rotmat, checks=false)
        @test abs(se2_from_so2_disp(r, disp).dc - dc) ≈ 0 atol=1e-16
    end

    r = so2_from_angle(0)
    for disp = ([1], [1, 2, 3], [1, 2, 3, 4])
        @test_throws DomainError se2_from_so2_disp(r, disp)
    end
end

@testset "one" begin
    for dc = good_dual_complexes
        @test one(se2_from_dual_complex(dc)).dc == one(DualComplex)
    end

    @test one(SE2).dc == one(DualComplex) 
end

@testset "dual_complex" begin
    for dc = good_dual_complexes
        @test dual_complex(SE2(dc)) == dc
    end
end

@testset "matrix" begin
    for mat = good_matrices
        @test matrix(se2_from_matrix(mat)) ≈ mat
    end
end

@testset "so2" begin
    for mat = good_matrices
        @test rotmat(so2(se2_from_matrix(mat))) ≈ mat[1:2, 1:2]
    end
end

@testset "disp" begin
    for ad = angles_disps
        d = ad.disp
        dc = angle_disp_to_dual_complex(ad.angle, ad.disp)
        @test disp(se2_from_dual_complex(dc)) ≈ d atol=1e-15
    end
end

@testset "multiplication" begin
    for ad = angles_disps
        θ, d = ad
        m = angle_disp_to_matrix(θ, d)
        dc = angle_disp_to_dual_complex(θ, d)
        @test matrix(one(SE2) * se2_from_dual_complex(dc)) ≈ m
        @test matrix(se2_from_dual_complex(dc) * one(SE2)) ≈ m
    end

    for n = 1:length(angles_disps) - 1
        ad1, ad2 = angles_disps[n], angles_disps[n + 1]
        θ1, d1 = ad1
        θ2, d2 = ad2
        dc1 = angle_disp_to_dual_complex(θ1, d1)
        dc2 = angle_disp_to_dual_complex(θ2, d2)
        m1 = angle_disp_to_matrix(θ1, d1)
        m2 = angle_disp_to_matrix(θ2, d2)
        l1 = se2_from_matrix(m1)
        l2 = se2_from_matrix(m2)
        @test matrix(l1 * l2) ≈ m1 * m2
        @test matrix(l2 * l1) ≈ m2 * m1
    end
end

@testset "inverse" begin
    for dc = good_dual_complexes
        l = se2_from_dual_complex(dc)
        @test inv(l).dc == conj(dc)
        @test matrix(inv(l) * l) ≈ LA.I
        @test matrix(l * inv(l)) ≈ LA.I
    end
end

@testset "exponential map" begin
    for ad = angles_disps
        θ, d = ad
        @test matrix(exp(SE2, [θ, d...])) ≈ angle_disp_to_matrix(θ, d)
    end
end

@testset "logarithm" begin
    function test_log(l, θ, d)
        v = log(l)
        θl = v[1]
        dl = v[2:3]
        @test θl <= π
        @test θl >= -π
        @test sin(θl) ≈ sin(θ)
        @test cos(θl) ≈ cos(θ)
        @test dl ≈ d atol=1e-15
    end

    for ad = angles_disps
        dc = angle_disp_to_dual_complex(ad.angle, ad.disp)
        l = se2_from_dual_complex(dc)
        test_log(l, ad.angle, ad.disp)
        @test dual_complex(exp(SE2, log(l))) ≈ dc
    end
end
