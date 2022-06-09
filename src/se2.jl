struct SE2{T<:DualComplex} <: LieGroup
    dc::T
end

# Constructors
function from_dual_complex(T::Type{<:SE2}, dc::DualComplex; checks::Bool=true)
    if checks && !(abs2(dc) ≈ 1)
        throw(
            DomainError(
                dc,
                "Dual Complex number must have a norm of 1, got $(abs(dc))."))
    end
    return T(dc)
end

function from_matrix(T::Type{<:SE2}, t::Matrix{<:Real}; checks::Bool=true)
    if checks
        if size(t) != (3, 3)
            throw(
                DomainError(
                    t,
                    "Matrix must be of size (3, 3), got $(size(t))."))
        elseif t[3, :] != [0, 0, 1]
            throw(
                DomainError(
                    t,
                    "Matrix last row must be [0, 0, 1], got $(t[3, :])."))
        end
    end
    r = from_rotmat(SO2, t[1:2, 1:2], checks=checks)
    return from_so2_disp(SE2, r, t[1:2, 3], checks=checks)
end

function from_so2_disp(
            T::Type{<:SE2},
            r::SO2,
            d::Vector{R};
            checks::Bool=true
        ) where {R<:Real}
    if checks && size(d) != (2,)
        throw(
            DomainError(
                d,
                "Displacement vector must have size (2,), got $(size(d))."))
    end
    θ = angle(r)
    chθ = cos(θ/2)
    shθ = sin(θ/2)
    dcr = DualComplex(chθ, shθ, shθ, shθ)
    dcd = DualComplex{R}(one(R), zero(R), -d[2], d[1])
    return T(dcd * dcr)
end

Base.one(q::SE2{DualComplex{T}}) where {T<:Number} = SE2{DualComplex{T}}(
    DualComplex(one(T)))
Base.one(T::Type{<:SE2}) = T(one(DualComplex))
Base.one(T::Type{SE2{DualComplex{R}}}) where {R<:Number} = T(DualComplex(
    one(R)))

# Selectors
dual_complex(q::SE2) = q.dc

function matrix(q::SE2)
    dc = q.dc
    r = DualComplex(dc.a, dc.b, dc.b, dc.b)
    t = dc * conj(r)
    return vcat(hcat(rotmat(so2(q)), [t.d; -t.c]), [0 0 1])
end

so2(q::SE2) = from_complex(SO2, complex(q.dc) ^ 2)

function disp(q::SE2)
    dc = q.dc
    r = DualComplex(dc.a, dc.b, dc.b, dc.b)
    t = dc * conj(r)
    return [t.d, -t.c]
end

# Operators
Base.:*(q::SE2, p::SE2) = SE2(q.dc * p.dc)
Base.inv(q::T) where {T<:SE2} = T(conj(q.dc))
# The Lie algebra is represented as a three-element vector, containing one
# element corresponding to the rotation angle (in the range [-π, π]), and two
# elements corresponding to the translation vector, in this order.
Base.exp(T::Type{<:SE2}, v::Vector{<:Real}) = from_so2_disp(
    SE2,
    from_angle(SO2, v[1]),
    v[2:3],
    checks=false)
Base.log(q::SE2) = [angle(so2(q)), disp(q)...]
