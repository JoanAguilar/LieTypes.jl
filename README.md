# LieTypes.jl

Julia types for Lie groups.

The main goal is for the provided types to allow the writing of Lie-group-agnostic code.

The interface consists of only 4 methods<sup>1</sup>:

 - `*(::T, ::T)::T where {T<:LieGroup}`: element multiplication.

 - `inv(::T)::T where {T<:LieGroup}`: element inverse.

 - `exp(::Type{T}, ::Array)::T where {T<:LieGroup}`: exponemential map.

 - `log(<:LieGroup)::Array`: logarithm.

Each concrete type provides their own constructors and selectors, as well as some encoding of the Lie algebra (always represented as an `Array`).

Currently 6 types are provided (`LieScalar`, `LieVector`, `SO2`, `SE2`, `SO3`, and `SE3`), all related to rigid-body transforms.


## Support Types

3 support types are provided:

 - `DualComplex`, which implements what is sometimes referred as "dual complex numbers" or "planar quaternions" (see, for instance, [Applications of dual quaternions to 2D geometry](https://en.wikipedia.org/wiki/Applications_of_dual_quaternions_to_2D_geometry)). This type is used as the underlying representation of `SE2`, and as a valid interface (via the `se2_from_dual_complex` constructor and the `dual_complex` selector).

 - `Quaternion`, which implements quaternions. This type is used as the underlying representation of `SO3`, and as a valid interface (via the `so3_from_quaternion` constructor and `quaternion` selector).

 - `Dual`, which implements dual numbers. The type `Dual{Quaternion{Float64}}` is used as the underlying representation of `SE3`, and as a valid interface (via the `se3_from_dual_quaternion` constructor and `dual_quaternion` selector).

<sup>1</sup> The provided type annotations are only valid for non-parametric subtypes of `LieGroup`.
