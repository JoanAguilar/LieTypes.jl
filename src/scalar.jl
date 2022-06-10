struct LieScalar{T<:Number} <: LieGroup
    s::T
end

# Constructors
lie_scalar_from_number(n::Number) = LieScalar(n)
Base.one(l::T) where {T<:LieScalar} = T(0)
Base.one(T::Type{<:LieScalar}) = T(0)

# Selectors
number(l::LieScalar) = l.s

# Operators
Base.:*(l::LieScalar, m::LieScalar) = LieScalar(l.s + m.s)
Base.inv(l::T) where {T<:LieScalar} = T(-l.s)
# The Lie algebra is represented as a scalar array, containing the value of the
# transformation. The choice of a scalar array rather than a scalar is to keep
# type consistency with other Lie algebras (all represented by arrays).
Base.exp(T::Type{<:LieScalar}, a::Array{<:Number, 0}) = T(a[])
Base.log(l::LieScalar) = fill(l.s)
