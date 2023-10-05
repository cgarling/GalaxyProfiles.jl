################################################################
### Default implementations for random sampling

"""
    sample2D_r!([rng::Random.AbstractRNG=Random.default_rng()], d::AbstractMassProfile, x::AbstractArray{<:Real})

Draw multiple samples of the radius following the surface density profile of `d` and write them to `x` in place. By default, falls back to [`quantile2D.(d, rand!(rng, x))`](@ref quantile2D). See also [`sample2D_r`](@ref).
"""
function sample2D_r!(rng::AbstractRNG, d::AbstractMassProfile, x::AbstractArray{<:Real})
    rand!(rng, x)
    x .= quantile2D.(d, x)
end
sample2D_r!(d::AbstractMassProfile, x::AbstractArray{<:Real}) = rand!(default_rng(), d, x)

# Non-Mutating forms
"""
    sample2D_r([rng::Random.AbstractRNG=Random.default_rng()], d::AbstractMassProfile [, dims...])

Draws samples of the radius following the surface density profile of `d` with shape `[dims...]`, which can either be a `Vararg{Int,N}` (e.g., `sample2D_r(rng, d, 2, 2)` to generate a 2x2 Matrix) or a `Dims{N}`, which is an `NTuple` of `N` `Int`s (e.g., `sample2D_r(rng, d, (2,2))`. By default, falls back to [`quantile2D.(d, rand(rng [, dims...]))`](@ref quantile2D). See also [`sample2D_r!`](@ref).
"""
sample2D_r(rng::AbstractRNG, d::AbstractMassProfile, dims::Vararg{Int, N}) where N =
    ( x = rand(rng, dims...); x .= quantile2D.(d, x) )
sample2D_r(d::AbstractMassProfile, dims::Vararg{Int, N}) where N =
    sample2D_r(default_rng(), d, dims...)
sample2D_r(rng::AbstractRNG, d::AbstractMassProfile, dims::Dims) = sample2D_r(rng, d, dims...)
sample2D_r(d::AbstractMassProfile, dims::Dims) = sample2D_r(default_rng(), d, dims)
sample2D_r(rng::AbstractRNG, d::AbstractMassProfile) = quantile2D(d, rand(rng))
sample2D_r(d::AbstractMassProfile) = sample2D_r(default_rng(), d)

#### Now do the 3D samplers
# Mutating forms

"""
    sample3D_r!([rng::Random.AbstractRNG=Random.default_rng()], d::AbstractDensity, x::AbstractArray{<:Real})

Draw multiple samples of the radius following the density profile of `d` and write them to `x` in place. By default, falls back to [`quantile3D.(d, rand!(rng, x))`](@ref quantile3D). See also [`sample3D_r`](@ref).
"""
sample3D_r!(rng::AbstractRNG, d::AbstractDensity, x::AbstractArray{<:Real}) = 
    ( rand!(rng, x); x .= quantile3D.(d, x) )
sample3D_r!(d::AbstractDensity, x::AbstractArray{<:Real}) = sample3D_r!(default_rng(), d, x)

# Non-mutating forms
"""
    sample3D_r([rng::Random.AbstractRNG=Random.default_rng()], d::AbstractDensity [, dims...])

Sample random points following the density profile of `d` with shape `[dims...]`, which can either be a `Vararg{Int,N}` (e.g., `sample3D_r(rng, d, 2, 2)` to generate a 2x2 Matrix) or a `Dims{N}`, which is an `NTuple` of `N` `Int`s (e.g., `sample3D_r(rng, d, (2,2))`. By default, falls back to [`quantile3D.(d, rand(rng [, dims...]))`](@ref quantile3D). See also [`sample3D_r!`](@ref).
"""
sample3D_r(rng::AbstractRNG, d::AbstractDensity, dims::Vararg{Int, N}) where N =
    ( x = rand(rng, dims...); x .= quantile3D.(d, x) )
sample3D_r(d::AbstractDensity, dims::Vararg{Int, N}) where N =
    sample3D_r(default_rng(), d, dims...)
sample3D_r(rng::AbstractRNG, d::AbstractDensity, dims::Dims) = sample3D_r(rng, d, dims...)
sample3D_r(d::AbstractDensity, dims::Dims) = sample3D_r(default_rng(), d, dims)
sample3D_r(rng::AbstractRNG, d::AbstractDensity) = quantile3D(d, rand(rng))
sample3D_r(d::AbstractDensity) = sample3D_r(default_rng(), d)
