function quaternion_to_rotmat(q)
    θ = 2 * acos(scalar(q))
    if θ ≈ 0 || θ ≈ 2 * π
        return [1 0 0; 0 1 0; 0 0 1]
    end

    v = LA.normalize(vector(q))
    skew(x) = [0 -x[3] x[2]; x[3] 0 -x[1]; -x[2] x[1] 0]

    return cos(θ) * LA.I + (1 - cos(θ)) * v * transpose(v) + sin(θ) * skew(v)
end

function quaternion_to_rotvec(q)
    θ = 2 * acos(scalar(q))
    if θ ≈ 0 || θ ≈ 2 * π
        return [0, 0, 0]
    else
        v = LA.normalize(vector(q))
        return θ * v
    end
end

function rotvec_to_axis_angle(v)
    θ = LA.norm(v)
    if θ ≈ 0
        return ([1, 0, 0], 0)
    else
        return (v / θ, θ)
    end
end

good_quaternions = (
    Quaternion(1, 0, 0, 0),
    Quaternion(0, 1, 0, 0),
    Quaternion(0, 0, 1, 0),
    Quaternion(0, 0, 0, 1),
    Quaternion(-1, 0, 0, 0),
    Quaternion(0, -1, 0, 0),
    Quaternion(0, 0, -1, 0),
    Quaternion(0, 0, 0, -1),
    Quaternion(√2/2, √2/(2*√3), √2/(2*√3), √2/(2*√3)),
    Quaternion(√2/2, -√2/(2*√3), -√2/(2*√3), -√2/(2*√3)),
    Quaternion(√2/2, √2/2, 0.0, 0.0),
    Quaternion(0.5, 0.0, √3/2, 0.0),
    Quaternion(-√3/2, 0.0, 0.0, 0.5))
bad_quaternions = (
    Quaternion(0, 0, 0, 0),
    Quaternion(2, 0, 0, 0),
    Quaternion(1, 1, 0, 0),
    Quaternion(1, 0, -1, 0),
    Quaternion(-1, 0, 0, 1),
    Quaternion(1, 2, 3, 4))
quaternions = (good_quaternions..., bad_quaternions...)

good_rotmats = map(quaternion_to_rotmat, good_quaternions)
bad_rotmats = (
    [1 1 1; 0 1 1; 0 0 1],
    zeros(3, 3),
    zeros(4, 4),
    ones(3, 3),
    -ones(3, 3))
rotmats = (good_rotmats..., bad_rotmats...)

good_rotvecs = map(quaternion_to_rotvec, good_quaternions)
bad_rotvecs = ([1, 2, 3, 4],)
rotvecs = (good_rotvecs..., bad_rotvecs...)

good_axes_angles = map(rotvec_to_axis_angle, good_rotvecs)
bad_axes_angles = (
    ([0, 0, 0], -2),
    ([2, 0, 0], -1),
    ([1, 1, 0], 0),
    ([-1, -1, -1], 1))
axes_angles = (good_axes_angles..., bad_axes_angles...)


@testset "constructor" begin
    for q = quaternions 
        @test SO3(q).q == q
    end
end

@testset "so3_from_quaternion" begin
    for q = quaternions
        @test so3_from_quaternion(q, checks=false).q == q
    end

    for q = good_quaternions
        @test so3_from_quaternion(q).q == q
    end

    for q = bad_quaternions
        @test_throws DomainError so3_from_quaternion(q)
    end
end

@testset "so3_from_rotmat" begin
    for r = rotmats
        @test_nowarn so3_from_rotmat(r, checks=false)
    end

    for n = 1:length(good_quaternions)
        q = good_quaternions[n]
        r = quaternion_to_rotmat(q)
        @test so3_from_rotmat(r).q ≈ q || so3_from_rotmat(r).q ≈ -q
    end

    for r = bad_rotmats
        @test_throws DomainError so3_from_rotmat(r)
    end
end

@testset "so3_from_axis_angle" begin
    for vθ = axes_angles
        v, θ = vθ
        @test_nowarn so3_from_axis_angle(v, θ, checks=false)
    end

    for n = 1:length(good_quaternions)
        q = good_quaternions[n]
        r = quaternion_to_rotvec(q)
        v, θ = rotvec_to_axis_angle(r)
        @test so3_from_axis_angle(v, θ).q ≈ q ||
            so3_from_axis_angle(v, θ).q ≈ -q
    end

    for vθ = bad_axes_angles
        v, θ = vθ
        @test_throws DomainError so3_from_axis_angle(v, θ)
    end
end

@testset "so3_from_rotvec" begin
    for n = 1:length(good_quaternions)
        q = good_quaternions[n]
        v = quaternion_to_rotvec(q)
        @test so3_from_rotvec(v).q ≈ q || so3_from_rotvec(v).q ≈ -q
    end

    for v = bad_rotvecs
        @test_throws DomainError so3_from_rotvec(v)
    end
end

@testset "one" begin
    for q = good_quaternions 
        @test one(so3_from_quaternion(q)).q == one(Quaternion)
    end

    @test one(SO3).q == one(Quaternion)
end

@testset "quaternion" begin
    for q = good_quaternions
        @test quaternion(SO3(q)) == q
    end
end

@testset "rotmat" begin
    for r = good_rotmats
        @test rotmat(so3_from_rotmat(r)) ≈ r
    end
end

@testset "axis" begin
    for vθ = good_axes_angles
        v, θ = vθ
        if !(θ % 2π ≈ 0)
            if abs(θ % 2π) >= π
                @test axis(so3_from_axis_angle(v, θ)) ≈ -v
            else
                @test axis(so3_from_axis_angle(v, θ)) ≈ v
            end
        end
    end
end

@testset "angle" begin
    for vθ = good_axes_angles
        v, θ = vθ
        @test angle(so3_from_axis_angle(v, θ)) >= 0
        @test angle(so3_from_axis_angle(v, θ)) <= π
        if abs(θ % 2π) > π
            @test angle(so3_from_axis_angle(v, θ)) == 2π - abs(θ % 2π)
        else
            @test angle(so3_from_axis_angle(v, θ)) == abs(θ % 2π)
        end
    end
end

@testset "rotvec" begin
    for q = good_quaternions
        l = so3_from_quaternion(q)
        if angle(l) > 0
            @test LA.normalize(rotvec(l)) ≈ axis(l)
        else
            @test rotvec(l) == [0, 0, 0]
        end
        @test LA.norm(rotvec(l)) >= 0
        @test LA.norm(rotvec(l)) <= π
        @test LA.norm(rotvec(l)) ≈ angle(l)
    end
end

@testset "multiplication" begin
    for q = good_quaternions
        @test quaternion(one(SO3) * so3_from_quaternion(q)) == q
        @test quaternion(so3_from_quaternion(q) * one(SO3)) == q
    end

    for n = 1:length(good_quaternions) - 1
        q1, q2 = good_quaternions[n], good_quaternions[n + 1]
        r1, r2 = quaternion_to_rotmat(q1), quaternion_to_rotmat(q2)
        l1, l2 = so3_from_quaternion(q1), so3_from_quaternion(q2)
        @test rotmat(l1 * l2) ≈ r1 * r2
    end
end

@testset "inverse" begin
    for q = good_quaternions
        r = quaternion_to_rotmat(q)
        l = so3_from_quaternion(q)
        @test rotmat(inv(l)) ≈ transpose(r)
        @test rotmat(inv(l) * l) ≈ LA.I
        @test rotmat(l * inv(l)) ≈ LA.I
    end
end

@testset "exponential map" begin
    for q = good_quaternions
        v = quaternion_to_rotvec(q)
        r = quaternion_to_rotmat(q)
        @test rotmat(exp(SO3, v)) ≈ r
    end
end

@testset "logarithm" begin
    for q = good_quaternions
        v = quaternion_to_rotvec(q)
        r = quaternion_to_rotmat(q)
        l = so3_from_quaternion(q)
        @test LA.norm(log(l)) <= π
        @test LA.norm(log(l)) ≈ angle(l)
        if angle(l) > 0
            @test LA.normalize(log(l)) ≈ axis(l)
        end
        @test rotmat(exp(SO3, log(l))) ≈ r
    end
end
