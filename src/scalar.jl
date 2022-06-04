struct LieScalar{T<:Number} <: LieGroup
    s::T
end

# Constructors
from_number(T::Type{<:LieScalar}, n::Number) = T(n)
Base.one(l::T) where {T<:LieScalar} = T(0)
Base.one(T::Type{<:LieScalar}) = T(0)

# Selectors
number(l::LieScalar) = l.s

# Operators
Base.:*(l::LieScalar, m::LieScalar) = LieScalar(l.s + m.s)
Base.inv(l::T) where {T<:LieScalar} = T(-l.s)
Base.exp(T::Type{<:LieScalar}, a::Array{<:Number, 0}) = T(a[])
Base.log(l::LieScalar) = fill(l.s)
