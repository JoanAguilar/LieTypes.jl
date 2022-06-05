struct SO2{T<:Complex} <: LieGroup
    c::T
end

# Contsructors
function from_complex(T::Type{<:SO2}, c::Complex; checks::Bool=true)
    if checks && !(abs2(c) ≈ 1)
        throw(
            DomainError(
                c,
                "Complex number must have a norm of 1, got $(abs(c))."))
    end
    return T(c)
end

from_angle(T::Type{<:SO2}, θ::Real) = T(cos(θ) + im * sin(θ))

function from_rotmat(T::Type{<:SO2}, r::Matrix{<:Real}; checks::Bool=true)
    if checks
        if size(r) != (2, 2)
            throw(
                DomainError(
                    r,
                    "Rotation matrix must be of size (2, 2), got $(size(r))."))
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
    return T(r[1, 1] + im * r[2, 1])
end

Base.one(q::T) where {T<:SO2} = T(Complex(1))
Base.one(T::Type{<:SO2}) = T(Complex(1))

# Selectors
Base.complex(q::SO2) = q.c
Base.angle(q::SO2) = angle(q.c)
rotmat(q::SO2) = [
    real(q.c) -imag(q.c);
    imag(q.c) real(q.c)]

# Operators
Base.:*(q::SO2, p::SO2) = SO2(q.c * p.c)
Base.inv(q::T) where {T<:SO2} = T(conj(q.c))
Base.exp(T::Type{<:SO2}, v::Array{<:Real, 0}) = from_angle(T, v[])
Base.log(q::SO2) = fill(angle(q.c))