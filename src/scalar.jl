struct LieScalar{T<:Number} <: LieGroup
    s::T
end

# Constructors
lie_scalar_from_number(n::T) where {T<:Number} = LieScalar{T}(n)
Base.one(l::LieScalar{T}) where {T<:Number} = LieScalar{T}(zero(T))
Base.one(T::Type{LieScalar{R}}) where {R<:Number}= T(zero(R))

# Selectors
number(l::LieScalar) = l.s

# Operators
Base.:*(l::LieScalar{T}, m::LieScalar{R}) where {T<:Number, R<:Number} =
    LieScalar{promote_type(T, R)}(l.s + m.s)
Base.inv(l::T) where {T<:LieScalar} = T(-l.s)
# The Lie algebra is represented as a scalar array, containing the value of the
# transformation. The choice of a scalar array rather than a scalar is to keep
# type consistency with other Lie algebras (all represented by arrays).
Base.exp(T::Type{LieScalar{R}}, a::Array{<:Number, 0}) where {R<:Number} =
    T(R(a[]))
Base.log(l::LieScalar) = fill(l.s)
