struct SO2{T<:Complex} <: LieGroup
    c::T
end

# Constructors
function so2_from_complex(c::Complex; checks::Bool=true)
    if checks && !(abs2(c) ≈ 1)
        throw(
            DomainError(
                c,
                "Complex number must have a norm of 1, got $(abs(c))."))
    end
    return SO2(c)
end

so2_from_angle(θ::Real) = SO2(cos(θ) + im * sin(θ))

function so2_from_rotmat(r::Matrix{<:Real}; checks::Bool=true)
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
    return SO2(r[1, 1] + im * r[2, 1])
end

Base.one(q::SO2{Complex{T}}) where {T<:Number} = SO2{Complex{T}}(
    Complex(one(T)))
Base.one(T::Type{<:SO2}) = T(one(Complex))
Base.one(T::Type{SO2{Complex{R}}}) where {R<:Number} = T(Complex(one(R)))

# Selectors
Base.complex(q::SO2) = q.c
Base.angle(q::SO2) = angle(q.c)
rotmat(q::SO2) = [
    real(q.c) -imag(q.c);
    imag(q.c) real(q.c)]

# Operators
Base.:*(q::SO2, p::SO2) = SO2(q.c * p.c)
Base.inv(q::T) where {T<:SO2} = T(conj(q.c))
# The Lie algebra is represented as a scalar array, containing the angle of the
# transformation in the range [-π, π].The choice of a scalar array rather than
# a scalar is to keep type consistency with other Lie algebras (all represented
# by arrays).
Base.exp(T::Type{<:SO2}, v::Array{<:Real, 0}) = so2_from_angle(v[])
Base.log(q::SO2) = fill(angle(q.c))
