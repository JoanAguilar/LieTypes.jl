function rotvec_disp_to_dual_quaternion(rotvec, disp)
    θ = LA.norm(rotvec)
    if θ > 0
        v = rotvec / θ
        θh = 0.5 * θ
        qr = Quaternion(cos(θh), sin(θh) * v...)
    else
        qr = one(Quaternion)
    end
    t = 0.5 * Quaternion(0, disp...) * qr
    return Dual(qr, t)
end

function rotvec_disp_to_matrix(rotvec, disp)
    return vcat(hcat(rotmat(so3_from_rotvec(rotvec)), disp), [0 0 0 1])
end

rotvecs_disps = (
    (rotvec=[0, 0, 0], disp=[0, 0, 0]),
    (rotvec=[π/6, 0, 0], disp=[0, 0, 0]),
    (rotvec=[0, 0.75π, 0], disp=[0, 0, 0]),
    (rotvec=[0, 0, -π/6], disp=[0, 0, 0]),
    (rotvec=[-0.75π, 0, 0], disp=[0, 0, 0]),
    (rotvec=[1, 1, 0], disp=[0, 0, 0]),
    (rotvec=[1, 1, 1], disp=[0, 0, 0]),
    (rotvec=[0, 0, 0], disp=[1, 0, 0]),
    (rotvec=[0, 0, 0], disp=[0, 1, 0]),
    (rotvec=[0, 0, 0], disp=[0, 1, 2]),
    (rotvec=[π/6, 0, 0], disp=[0, -1, 2]),
    (rotvec=[0, -π/3, 0], disp=[3, -4, 0]),
    (rotvec=[1, 1, 1], disp=[2, 3, -4]))

rotvecs_disps_good_dual_quaternions = map(
    rd -> rotvec_disp_to_dual_quaternion(rd.rotvec, rd.disp),
    rotvecs_disps)
other_good_dual_quaternions = (
   Dual(Quaternion(1), Quaternion(0, 2, 3, 4)),
   Dual(Quaternion(0, 1, 0, 0), Quaternion(0, 0, 0, 0)),
   Dual(Quaternion(0, 0, 1, 0), Quaternion(0, 0, 0, 0)),
   Dual(Quaternion(0, 0, 0, 1), Quaternion(0, 0, 0, 0)))
good_dual_quaternions = (
    rotvecs_disps_good_dual_quaternions...,
    other_good_dual_quaternions...)
bad_dual_quaternions = (
    zero(Dual{Quaternion}),
    Dual{Quaternion}(Quaternion(2), zero(Quaternion)),
    Dual{Quaternion}(Quaternion(0, -2, 0, 0), Quaternion(1, 2, 3, 4)))
dual_quaternions = (good_dual_quaternions..., bad_dual_quaternions...)

good_matrices = map(
    rd -> rotvec_disp_to_matrix(rd.rotvec, rd.disp),
    rotvecs_disps)
bad_matrices = (
    [1 1 1 2;
     0 1 1 3;
     0 0 1 4;
     0 0 0 1],
    [1 0 0 2;
     0 1 0 3;
     0 0 1 4;
     5 6 7 8],
    zeros(4, 4),
    ones(4, 4),
    zeros(5, 5))
matrices = (good_matrices..., bad_matrices...)

@testset "constructor" begin
    for dq = dual_quaternions
        @test SE3(dq).dq == dq
    end
end

@testset "se3_from_dual_quaternion" begin
    for dq = dual_quaternions
        @test se3_from_dual_quaternion(dq, checks=false).dq == dq
    end

    for dq = good_dual_quaternions
        @test se3_from_dual_quaternion(dq).dq == dq
    end

    for dq = bad_dual_quaternions
        @test_throws DomainError se3_from_dual_quaternion(dq)
    end
end

@testset "se3_from_matrix" begin
    for mat in matrices
        @test_nowarn se3_from_matrix(mat, checks=false)
    end

    for rd in rotvecs_disps
        v, d = rd
        dq = rotvec_disp_to_dual_quaternion(v, d)
        mat = rotvec_disp_to_matrix(v, d)
        @test se3_from_matrix(mat).dq ≈ dq || se3_from_matrix(mat).dq ≈ -dq
    end

    for mat in bad_matrices
        @test_throws DomainError se3_from_matrix(mat)
    end
end

@testset "se3_from_so3_disp" begin
    for mat in matrices
        rotmat = mat[1:3, 1:3]
        r = so3_from_rotmat(rotmat, checks=false)
        disp = mat[1:3, 4]
        @test_nowarn se3_from_so3_disp(r, disp, checks=false)
    end

    for rd in rotvecs_disps
        v, d = rd
        dq = rotvec_disp_to_dual_quaternion(v, d)
        mat = rotvec_disp_to_matrix(v, d)
        rotmat = mat[1:3, 1:3]
        r = so3_from_rotmat(rotmat, checks=false)
        @test (
            se3_from_so3_disp(r, d).dq ≈ dq ||
            se3_from_so3_disp(r, d).dq ≈ -dq)
    end

    r = so3_from_rotvec([0, 0, 0])
    for d = ([1], [1, 2], [1, 2, 3, 4])
        @test_throws DomainError se3_from_so3_disp(r, d)
    end
end

@testset "one" begin
    for dq = good_dual_quaternions
        @test one(se3_from_dual_quaternion(dq)).dq == one(Dual{Quaternion})
    end

    @test one(SE3).dq == one(Dual{Quaternion}) 
    @test one(SE3{Dual{Quaternion}}).dq == one(Dual{Quaternion})
    for T = (Int64, Float64)
        @test one(SE3{Dual{Quaternion{T}}}).dq == one(Dual{Quaternion})
    end
end

@testset "dual_quaternion" begin
    for dq = good_dual_quaternions
        @test dual_quaternion(SE3(dq)) == dq
    end
end

@testset "matrix" begin
    for mat = good_matrices
        @test matrix(se3_from_matrix(mat)) ≈ mat
    end
end

@testset "so3" begin
    for mat = good_matrices
        @test rotmat(so3(se3_from_matrix(mat))) ≈ mat[1:3, 1:3]
    end
end

@testset "disp" begin
    for rd = rotvecs_disps
        d = rd.disp
        dq = rotvec_disp_to_dual_quaternion(rd.rotvec, rd.disp)
        @test disp(se3_from_dual_quaternion(dq)) ≈ d atol=1e-14
    end
end

@testset "multiplication" begin
    for rd = rotvecs_disps
        v, d = rd
        m = rotvec_disp_to_matrix(v, d)
        dq = rotvec_disp_to_dual_quaternion(v, d)
        @test matrix(one(SE3) * se3_from_dual_quaternion(dq)) ≈ m
        @test matrix(se3_from_dual_quaternion(dq) * one(SE3)) ≈ m
    end

    for n = 1:length(rotvecs_disps) - 1
        rd1, rd2 = rotvecs_disps[n], rotvecs_disps[n + 1]
        v1, d1 = rd1
        v2, d2 = rd2
        dq1 = rotvec_disp_to_dual_quaternion(v1, d1)
        dq2 = rotvec_disp_to_dual_quaternion(v2, d2)
        m1 = rotvec_disp_to_matrix(v1, d1)
        m2 = rotvec_disp_to_matrix(v2, d2)
        l1 = se3_from_matrix(m1)
        l2 = se3_from_matrix(m2)
        @test matrix(l1 * l2) ≈ m1 * m2
        @test matrix(l2 * l1) ≈ m2 * m1
    end
end

@testset "inverse" begin
    for dq = good_dual_quaternions
        l = se3_from_dual_quaternion(dq)
        @test inv(l).dq == conj(dq)
        @test matrix(inv(l) * l) ≈ LA.I
        @test matrix(l * inv(l)) ≈ LA.I
    end
end

@testset "exponential map" begin
    for rd = rotvecs_disps
        v, d = rd
        @test matrix(exp(SE3, [v..., d...])) ≈ rotvec_disp_to_matrix(v, d)
    end
end

@testset "logarithm" begin
    function test_log(l, rv, d)
        v = log(l)
        rvl = v[1:3]
        dl = v[4:6]
        @test LA.norm(rvl) <= π
        @test isapprox(rvl / LA.norm(rvl), rv / LA.norm(rv), nans=true)
        @test cos(LA.norm(rvl)) ≈ cos(LA.norm(rv))
        @test dl ≈ d atol=1e-14
    end

    for rd = rotvecs_disps
        dq = rotvec_disp_to_dual_quaternion(rd.rotvec, rd.disp)
        l = se3_from_dual_quaternion(dq)
        test_log(l, rd.rotvec, rd.disp)
        @test dual_quaternion(exp(SE3, log(l))) ≈ dq
    end
end
