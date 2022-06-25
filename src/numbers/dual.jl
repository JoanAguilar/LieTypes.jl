struct Dual{T<:Number} <:Number
    r::T
    d::T
end

function Base.promote_rule(
            ::Type{<:Dual{S}},
            ::Type{<:Dual{U}}
        ) where {S<:Number, U<:Number}
    return Dual{promote_type(S, U)}
end

function Base.promote_rule(
            T::Type{<:Real},
            ::Type{<:Dual{S}}
        ) where {S<:Number}
    return Dual{promote_type(T, S)}
end

# Constructors
Dual(r::T; d::T=zero(T)) where {T<:Number} = Dual{T}(r, d)
Dual(r::Number, d::Number) = Dual(promote(r, d)...)
Base.zero(T::Type{Dual{R}}) where{R<:Number} = T(zero(R), zero(R))
Base.zero(T::Type{Dual}) = T(0)
Base.zero(n::T) where{T<:Dual} = zero(T)
Base.one(T::Type{Dual{R}}) where{R<:Number} = T(one(R), zero(R))
Base.one(T::Type{Dual}) = T(1)
Base.one(n::T) where{T<:Dual} = one(T)
du = Dual(0, 1)
ϵ = du

# Selectors
Base.real(n::Dual) = n.r
truereal(n::Dual) = real(n.r)
redu(n::Dual) = (n.r, n.d)
dual(n::Dual) = n.d
Base.isreal(n::Dual) = n.d == zero(n.d) && isreal(n.r)

# Operators
Base.abs2(n::Dual) = n * conj(n)
Base.conj(n::T) where {T<:Dual} = T(conj(n.r), conj(n.d))
dualconj(n::T) where {T<:Dual} = T(n.r, -n.d)

function Base.inv(n::Dual{R}) where {R<:Number}
    inv_real = inv(real(n))
    return Dual(inv_real, -inv_real * dual(n) * inv_real)
end

Base.:+(m::Dual, n::Dual) = Dual(m.r + n.r, m.d + n.d)
Base.:+(m::Number, n::Dual) = +(Dual(m), n)
Base.:+(m::Dual, n::Number) = +(n, m)
Base.:-(n::T) where {T<:Dual} = T(-n.r, -n.d)
Base.:-(m::Dual, n::Dual) = +(m, -n)
Base.:-(m::Number, n::Dual) = -(Dual(m), n)
Base.:-(m::Dual, n::Number) = -(m, Dual(n))
Base.:*(m::Dual, n::Dual) = Dual(m.r * n.r, m.r * n.d + m.d * n.r)
Base.:*(m::Number, n::Dual) = *(Dual(m) * n)
Base.:*(m::Dual, n::Number) = *(m * Dual(n))
Base.:/(m::Dual, n::Dual) = *(m, inv(n))
Base.:/(m::Number, n::Dual) = /(Dual(m), n)
Base.:/(m::Dual, n::Number) = /(m, Dual(n))
Base.:\(m::Dual, n::Dual) = *(inv(m), n)
Base.:\(m::Number, n::Dual) = \(Dual(m), n)
Base.:\(m::Dual, n::Number) = \(m, Dual(n))
Base.:(==)(m::Dual, n::Dual) = m.r == n.r && m.d == n.d
Base.:(==)(m::Dual, n::Number) = (==)(m, Dual(n))
Base.:(==)(m::Number, n::Dual) = (==)(n, m)

function dualrtoldefault(m::Dual, n::Dual, atol::Real)
    rtol = max(
        Base.rtoldefault(typeof(truereal(m))),
        Base.rtoldefault(typeof(truereal(n))))
    return atol > 0 ? zero(rtol) : rtol
end

dualrtoldefault(m::Number, n::Dual, atol::Real) =
    dualrtoldefault(Dual(m), n, atol)
dualrtoldefault(m::Dual, n::Number, atol::Real) =
    dualrtoldefault(m, Dual(n), atol)

function Base.isapprox(
        m::Dual,
        n::Dual;
        atol::Real=0,
        rtol::Real=dualrtoldefault(m, n, atol),
        nans::Bool=false,
        norm::Function=abs)
    return isapprox(m.r, n.r, atol=atol, rtol=rtol, nans=nans, norm=norm) &&
        isapprox(m.d, n.d, atol=atol, rtol=rtol, nans=nans, norm=norm)
end

Base.isapprox(
    m::Number,
    n::Dual;
    atol::Real=0,
    rtol::Real=dualrtoldefault(m, n, atol),
    nans::Bool=false,
    norm::Function=abs
) = isapprox(Dual(m), n, atol=atol, rtol=rtol, nans=nans, norm=norm)
Base.isapprox(
    m::Dual,
    n::Number;
    atol::Real=0,
    rtol::Real=dualrtoldefault(m, n, atol),
    nans::Bool=false,
    norm::Function=abs
) = isapprox(m, Dual(n), atol=atol, rtol=rtol, nans=nans, norm=norm)
