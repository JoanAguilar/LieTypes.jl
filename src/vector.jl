struct LieVector{T<:Array} <: LieGroup
    a::T
end

# Constructors
from_array(T::Type{<:LieVector}, a::Array) = T(a)
Base.one(l::T) where {T<:LieVector} = T(zero(l.a))

# Selectors
array(l::LieVector) = l.a

# Operators
Base.:*(l::LieVector, m::LieVector) = LieVector(l.a + m.a)
Base.inv(l::T) where {T<:LieVector} = T(-l.a)
Base.exp(T::Type{<:LieVector}, a::Array) = T(a)
Base.log(l::LieVector) = l.a
