struct DualComplex{T<:Number} <: Number
    a::T
    b::T
    c::T
    d::T
end

# Constructors
function DualComplex(
            a::T;
            b::T=zero(T),
            c::T=zero(T),
            d::T=zero(T)
        ) where {T<:Number}
    return DualComplex{T}(a, b, c, d)
end

DualComplex(c::Complex{T}) where {T<:Number} = DualComplex{T}(
    real(c),
    imag(c),
    zero(T),
    zero(T))
Base.zero(T::Type{DualComplex{R}}) where {R<:Number} = T(
    zero(R),
    zero(R),
    zero(R),
    zero(R))
Base.zero(T::Type{DualComplex}) = T(0)
Base.zero(dc::T) where {T<:DualComplex} = zero(T)
Base.one(T::Type{DualComplex{R}}) where {R<:Number} = T(
    one(R),
    zero(R),
    zero(R),
    zero(R))
Base.one(T::Type{DualComplex}) = T(1)
Base.one(dc::T) where {T<:DualComplex} = one(T)

# Selectors
scalar(dc::DualComplex) = dc.a
Base.real(dc::DualComplex) = scalar(dc)
Base.complex(dc::DualComplex{T}) where {T<:Number} = Complex{T}(dc.a, dc.b)
Base.imag(dc::DualComplex) = dc.b
vector(dc::DualComplex) = [dc.b, dc.c, dc.d]
dual(dc::DualComplex) = [dc.c, dc.d]

# Operators
Base.abs2(dc::DualComplex) = abs2(dc.a) + abs2(dc.b)
Base.abs(dc::DualComplex) = sqrt(abs2(dc))
Base.conj(dc::T) where {T<:DualComplex} = T(dc.a, -dc.b, -dc.c, -dc.d)
Base.inv(dc::DualComplex) = conj(dc) / abs2(dc)
Base.isreal(dc::DualComplex) = dc.b == 0 && dc.c == 0 && dc.d == 0
iscomplex(dc::DualComplex) = dc.c == 0 && dc.d == 0
Base.:+(u::DualComplex, v::DualComplex) = DualComplex(
    u.a + v.a,
    u.b + v.b,
    u.c + v.c,
    u.d + v.d)
Base.:+(u::Real, v::DualComplex) = DualComplex(u) + v
Base.:+(u::DualComplex, v::Real) = v + u
Base.:+(u::Complex, v::DualComplex) = DualComplex(u) + v 
Base.:+(u::DualComplex, v::Complex) = v + u
Base.:-(u::T) where{T<:DualComplex} = T(-u.a, -u.b, -u.c, -u.d)
Base.:-(u::DualComplex, v::DualComplex) = +(u, -v)
Base.:-(u::Real, v::DualComplex) = -(DualComplex(u), v)
Base.:-(u::DualComplex, v::Real) = -(u, DualComplex(v))
Base.:-(u::Complex, v::DualComplex) = -(DualComplex(u), v)
Base.:-(u::DualComplex, v::Complex) = -(u, DualComplex(v))
Base.:*(u::DualComplex, v::DualComplex) = DualComplex(
    u.a * v.a - u.b * v.b,
    u.a * v.b + u.b * v.a,
    u.a * v.c - u.b * v.d + u.c * v.a + u.d * v.b,
    u.a * v.d + u.b * v.c - u.c * v.b + u.d * v.a)
Base.:*(u::Real, v::DualComplex) = *(DualComplex(u) * v)
Base.:*(u::DualComplex, v::Real) = *(u * DualComplex(v))
Base.:*(u::Complex, v::DualComplex) = *(DualComplex(u) * v)
Base.:*(u::DualComplex, v::Complex) = *(u * DualComplex(v))
Base.:/(u::DualComplex, v::DualComplex) = *(u, inv(v))
Base.:/(u::Real, v::DualComplex) = /(DualComplex(u), v)
Base.:/(u::DualComplex, v::Real) = DualComplex(
    u.a / v,
    u.b / v,
    u.c / v,
    u.d / v)
Base.:/(u::Complex, v::DualComplex) = /(DualComplex(u), v)
Base.:/(u::DualComplex, v::Complex) = /(u, DualComplex(v))
Base.:\(u::DualComplex, v::DualComplex) = /(v, u)
Base.:\(u::Real, v::DualComplex) = /(v, u)
Base.:\(u::DualComplex, v::Real) = /(v, u)
Base.:\(u::Complex, v::DualComplex) = /(v, u)
Base.:\(u::DualComplex, v::Complex) = /(v, u)
Base.:(==)(u::DualComplex, v::DualComplex) =
    u.a == v.a && u.b == v.b && u.c == v.c && u.d == v.d
# Using `abs` as `norm`, and ignoring `p`, matches the behavior of `Complex`.
LA.norm(u::DualComplex, p::Real=2) = abs(u)
LA.normalize(u::DualComplex) = u / abs(u)
