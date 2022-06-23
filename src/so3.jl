struct SO3{T<:Quaternion} <: LieGroup
    q::T
end

# Constructors
function so3_from_quaternion(q::Quaternion; checks::Bool=true)
    if checks && !(abs2(q) ≈ 1)
        throw(
            DomainError(q, "Quaternion must have a norm of 1, got $(abs(q))."))
    end
    return SO3(q)
end

function so3_from_rotmat(r::Matrix{<:Real}; checks::Bool=true)
    if checks
        if size(r) != (3, 3)
            throw(
                DomainError(
                    r,
                    "Rotation matrix must be of size (3, 3), got $(size(r))."))
        elseif !(LA.det(r) ≈ 1)
            throw(
                DomainError(
                    r,
                    "Rotation matrix determinant must be 1, got " *
                    "$(LA.det(r))."))
        elseif !(transpose(r) * r ≈ LA.I && r * transpose(r) ≈ LA.I)
            throw(
                DomainError(
                    r,
                    "The transpose of the rotation matrix must also be its " *
                    "inverse."))
        end
    end

    # Adapted from:
    # Converting a Rotation Matrix to a Quaternion
    # by Mike Day, Insomniac Games
    if r[3, 3] < 0
        if r[1, 1] > r[2, 2]
            t = 0.5 * sqrt(1 + r[1, 1] - r[2, 2] - r[3, 3])
            q = Quaternion(
                (r[3, 2] - r[2, 3]) / (4 * t),
                t,
                (r[2, 1] + r[1, 2]) / (4 * t),
                (r[1, 3] + r[3, 1]) / (4 * t))
        else
            t = 0.5 * sqrt(1 - r[1, 1] + r[2, 2] - r[3, 3])
            q = Quaternion(
                (r[1, 3] - r[3, 1]) / (4 * t),
                (r[2, 1] + r[1, 2]) / (4 * t),
                t,
                (r[3, 2] + r[2, 3]) / (4 * t))
        end
    else
        if r[1, 1] < -r[2, 2]
            t = 0.5 * sqrt(1 - r[1, 1] - r[2, 2] + r[3, 3])
            q = Quaternion(
                (r[2, 1] - r[1, 2]) / (4 * t),
                (r[1, 3] + r[3, 1]) / (4 * t),
                (r[3, 2] + r[2, 3]) / (4 * t),
                t)
        else
            t = 0.5 * sqrt(1 + r[1, 1] + r[2, 2] + r[3, 3])
            q = Quaternion(
                t,
                (r[3, 2] - r[2, 3]) / (4 * t),
                (r[1, 3] - r[3, 1]) / (4 * t),
                (r[2, 1] - r[1, 2]) / (4 * t))
        end
    end
    return SO3(q)
end

function so3_from_axis_angle(v::Vector{<:Real}, θ::Real; checks::Bool=true)
    if checks
        if size(v) != (3,)
            throw(
                DomainError(
                    v,
                    "Axis vector must have size (3,), got $(size(v))."))
        elseif !(LA.norm(v) ≈ 1)
            throw(
                DomainError(
                    v,
                    "Axis vector must have norm 1, got $(LA.norm(v))."))
        end
    end

    chθ = cos(0.5 * θ)
    shθ = sin(0.5 * θ)
    return SO3(Quaternion(chθ, shθ * v...))
end

function so3_from_rotvec(v::Vector{<:Real}; checks::Bool=true)
    if checks && size(v) != (3,)
        throw(
            DomainError(
                v,
                "Rotation vector must have size (3,), got $(size(v))."))
    end

    θ = LA.norm(v)
    if θ ≈ 0
        return one(SO3)
    else
        return so3_from_axis_angle(v / θ, θ, checks=false)
    end
end

Base.one(q::SO3{Quaternion{T}}) where {T<:Number} = SO3{Quaternion{T}}(
    Quaternion(one(T)))
Base.one(T::Type{SO3}) = T(one(Quaternion))
Base.one(T::Type{SO3{Quaternion}}) = T(one(Quaternion))
Base.one(T::Type{SO3{Quaternion{R}}}) where {R<:Number} = T(Quaternion(one(R)))

# Selectors
quaternion(q::SO3) = q.q

function rotmat(q::SO3)
    s, i, j, k = sijk(q.q)
    return [
        -1 + 2 * (i^2 + s^2) 2 * (i * j - k * s)   2 * (i * k + j * s);
        2 * (i * j + k * s)  -1 + 2 * (j^2 + s ^2) 2 * (j * k - i * s);
        2 * (i * k - j * s)  2 * (j * k + i * s)   -1 + 2 * (k^2 + s^2)]
end

function axis(q::SO3)
    θ = 2 * acos(scalar(q.q))
    if θ ≈ 0 || θ ≈ 2 * π
        return [0, 0, 0]
    end

    v = vector(q.q) / LA.norm(vector(q.q))

    if θ > π
        # Account for angle wraparound.
        return -v 
    else
        return v
    end
end


function Base.angle(q::SO3)
    θ = 2 * acos(scalar(q.q))

    # Account for angle wraparound.
    if θ > π
        θ = 2 * π - θ
    end

    return θ
end

function rotvec(q::SO3)
    θ = angle(q)
    if θ ≈ 0
        return [0, 0, 0]
    else
        return θ * axis(q)
    end
end

# Operators
Base.:*(q::SO3, p::SO3) = SO3(q.q * p.q)
Base.inv(q::T) where {T<:SO3} = T(conj(q.q))
# The Lie algebra is represented as a three-element vector, containing the
# rotation vector of the transformation.
Base.exp(T::Type{<:SO3}, v::Vector{<:Real}) = so3_from_rotvec(v, checks=false)
Base.log(q::SO3) = rotvec(q)
