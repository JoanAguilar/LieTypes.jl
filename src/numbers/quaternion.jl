struct Quaternion{T<:Number} <: Number
    s::T
    i::T
    j::T
    k::T
end

# Constructors
function Quaternion(
            s::T;
            i::T=zero(T),
            j::T=zero(T),
            k::T=zero(T)
        ) where {T<:Number}
    return Quaternion{T}(s, i, j, k)
end

Quaternion(c::Complex{T}) where {T<:Number} = Quaternion{T}(
    real(c),
    imag(c),
    zero(T),
    zero(T))
Base.zero(T::Type{Quaternion{R}}) where {R<:Number} = T(
    zero(R),
    zero(R),
    zero(R),
    zero(R))
Base.zero(T::Type{Quaternion}) = T(0)
Base.zero(q::T) where {T<:Quaternion} = zero(T)
Base.one(T::Type{Quaternion{R}}) where {R<:Number} = T(
    one(R),
    zero(R),
    zero(R),
    zero(R))
Base.one(T::Type{Quaternion}) = T(1)
Base.one(q::T) where {T<:Quaternion} = one(T)

# Selectors
scalar(q::Quaternion) = q.s
Base.real(q::Quaternion) = scalar(q)
Base.complex(q::Quaternion{T}) where {T<:Number} = Complex{T}(q.s, q.i)
sijk(q::Quaternion) = (q.s, q.i, q.j, q.k)
Base.imag(q::Quaternion) = q.i
vector(q::Quaternion) = [q.i, q.j, q.k]
Base.isreal(q::Quaternion) = q.i == 0 && q.j == 0 && q.k == 0
iscomplex(q::Quaternion) = q.j == 0 && q.k == 0

# Operators
Base.abs2(q::Quaternion) = abs2(q.s) + abs2(q.i) + abs2(q.j) + abs2(q.k)
Base.abs(q::Quaternion) = sqrt(abs2(q))
Base.conj(q::T) where {T<:Quaternion} = T(q.s, -q.i, -q.j, -q.k)
Base.inv(q::Quaternion) = conj(q) / abs2(q)
Base.:+(p::Quaternion, q::Quaternion) = Quaternion(
    p.s + q.s,
    p.i + q.i,
    p.j + q.j,
    p.k + q.k)
Base.:+(u::Real, q::Quaternion) = Quaternion(u) + q
Base.:+(p::Quaternion, v::Real) = v + p
Base.:+(u::Complex, q::Quaternion) = Quaternion(u) + q 
Base.:+(p::Quaternion, v::Complex) = v + p
Base.:-(q::T) where{T<:Quaternion} = T(-q.s, -q.i, -q.j, -q.k)
Base.:-(p::Quaternion, q::Quaternion) = +(p, -q)
Base.:-(u::Real, q::Quaternion) = -(Quaternion(u), q)
Base.:-(p::Quaternion, v::Real) = -(p, Quaternion(v))
Base.:-(u::Complex, q::Quaternion) = -(Quaternion(u), q)
Base.:-(p::Quaternion, v::Complex) = -(p, Quaternion(v))
Base.:*(p::Quaternion, q::Quaternion) =
    Quaternion(
        p.s * q.s - p.i * q.i - p.j * q.j - p.k * q.k,
        p.s * q.i + p.i * q.s + p.j * q.k - p.k * q.j,
        p.s * q.j - p.i * q.k + p.j * q.s + p.k * q.i,
        p.s * q.k + p.i * q.j - p.j * q.i + p.k * q.s)
Base.:*(u::Real, q::Quaternion) = *(Quaternion(u) * q)
Base.:*(p::Quaternion, v::Real) = *(p * Quaternion(v))
Base.:*(u::Complex, q::Quaternion) = *(Quaternion(u) * q)
Base.:*(p::Quaternion, v::Complex) = *(p * Quaternion(v))
Base.:/(p::Quaternion, q::Quaternion) = *(p, inv(q))
Base.:/(u::Real, q::Quaternion) = /(Quaternion(u), q)
Base.:/(p::Quaternion, v::Real) = Quaternion(
    p.s / v,
    p.i / v,
    p.j / v,
    p.k / v)
Base.:/(u::Complex, q::Quaternion) = /(Quaternion(u), q)
Base.:/(p::Quaternion, v::Complex) = /(p, Quaternion(v))
Base.:\(p::Quaternion, q::Quaternion) = /(q, p)
Base.:\(u::Real, q::Quaternion) = /(q, u)
Base.:\(p::Quaternion, v::Real) = /(v, p)
Base.:\(u::Complex, q::Quaternion) = /(q, u)
Base.:\(p::Quaternion, v::Complex) = /(v, p)
Base.:(==)(p::Quaternion, q::Quaternion) =
    q.s == p.s && q.i == p.i && q.j == p.j && q.k == p.k
# Using `abs` as `norm`, and ignoring `p`, matches the behavior of `Complex`.
LA.norm(q::Quaternion, p::Real=2) = abs(q)
LA.normalize(q::Quaternion) = q / abs(q)
