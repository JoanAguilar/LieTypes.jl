struct SE3{T<:Dual{<:Quaternion}} <: LieGroup
    dq::T
end

# Constructors
function se3_from_dual_quaternion(dq::Dual{<:Quaternion}; checks::Bool=true)
    # `abs2(::Dual{Quaternion})` returns a `Dual{Quaternion}`, and
    # `isapprox(::Dual{Quaternion}, ::Integer)` ends up comparing elements with
    # 0 magnitude, which means we need to raise `atol` in order for `isapprox`
    # to be more permissive than `==`. As `dq` is expected to have norm 1 in
    # the well-behaved cases, we use the default `rtol` value as `atol` to
    # circumvent the problem.
    if checks && !(isapprox(abs2(dq), 1, atol=dualrtoldefault(dq, 1, 0)))
        throw(
            DomainError(
                dq,
                "Dual Quaternion number must have a norm of 1, got a " *
                "squared norm of $(abs2(dq))."))
    end
    return SE3(dq)
end

function se3_from_matrix(t::Matrix{<:Real}; checks::Bool=true)
    if checks
        if size(t) != (4, 4)
            throw(
                DomainError(
                    t,
                    "Matrix must be of size (4, 4), got $(size(t))."))
        elseif t[4, :] != [0, 0, 0, 1]
            throw(
                DomainError(
                    t,
                    "Matrix last row must be [0, 0, 0, 1], got $(t[4, :])."))
        end
    end
    r = so3_from_rotmat(t[1:3, 1:3], checks=checks)
    return se3_from_so3_disp(r, t[1:3, 4], checks=checks)
end

function se3_from_so3_disp(
            r::SO3,
            d::Vector{R};
            checks::Bool=true
        ) where {R<:Real}
    if checks && size(d) != (3,)
        throw(
            DomainError(
                d,
                "Displacement vector must have size (3,), got $(size(d))."))
    end
    qr = quaternion(r)
    t = Quaternion(zero(R), d...)
    return SE3(Dual(qr, 0.5 * t * qr))
end

Base.one(q::SE3{Dual{Quaternion{T}}}) where {T<:Real} =
    SE3{Dual{Quaternion{T}}}(one(Dual{Quaternion{T}}))
Base.one(T::Type{SE3}) = T(one(Dual{Quaternion}))
Base.one(T::Type{SE3{Dual{Quaternion}}}) = T(one(Dual{Quaternion}))
Base.one(T::Type{SE3{Dual{Quaternion{R}}}}) where {R<:Real} =
    T(one(Dual{Quaternion{R}}))

# Selectors
dual_quaternion(q::SE3) = q.dq

function matrix(q::SE3)
    dq = q.dq
    qr = real(dq)
    t = vector(2 * dual(dq) * conj(qr))
    return vcat(
        hcat(rotmat(so3_from_quaternion(qr, checks=false)), t),
        [0 0 0 1])
end

so3(q::SE3) = so3_from_quaternion(real(q.dq))
disp(q::SE3) = vector(2 * dual(q.dq) * conj(real(q.dq)))

# Operators
Base.:*(q::SE3, p::SE3) = SE3(q.dq * p.dq)
Base.inv(q::T) where {T<:SE3} = T(conj(q.dq))
# The Lie algebra is represented as a six-element vector, containing three 
# elements corresponding to the rotation vector, and three elements
# corresponding to the translation vector, in this order.
Base.exp(T::Type{<:SE3}, v::Vector{<:Real}) = se3_from_so3_disp(
    so3_from_rotvec(v[1:3]),
    v[4:6],
    checks=false)
Base.log(q::SE3) = [rotvec(so3(q)); disp(q)...]
