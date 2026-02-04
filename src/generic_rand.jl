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

Draws samples of the radius following the surface density profile of `d` with shape `[dims...]`, which can either be a `Vararg{Integer, N}` (e.g., `sample2D_r(rng, d, 2, 2)` to generate a 2x2 Matrix) or a `Dims{N}`, which is an `NTuple` of `N` `Int`s (e.g., `sample2D_r(rng, d, (2,2))`. By default, falls back to [`quantile2D.(d, rand(rng [, dims...]))`](@ref quantile2D). See also [`sample2D_r!`](@ref).
"""
function sample2D_r(rng::AbstractRNG, d::AbstractMassProfile, dims::Integer...)
    return [sample2D_r(rng, d) for _ in CartesianIndices(dims)]
end
sample2D_r(rng::AbstractRNG, d::AbstractMassProfile, dims::Dims) = sample2D_r(rng, d, dims...)
sample2D_r(d::AbstractMassProfile, dims::Dims) = sample2D_r(default_rng(), d, dims)
sample2D_r(rng::AbstractRNG, d::AbstractMassProfile) = quantile2D(d, rand(rng))
sample2D_r(d::AbstractMassProfile) = sample2D_r(default_rng(), d)

"""
    sample2D(d::AbstractMassProfile, center::SkyCoords.AbstractSkyCoords, distance; 
             dims..., rng=Random.default_rng())

Sample 2D sky coordinates from the mass profile `d`.

# Arguments
- `d::AbstractMassProfile`: The mass profile to sample from
- `center::SkyCoords.AbstractSkyCoords`: The center position of the galaxy
- `distance`: Distance to the galaxy in the same physical units as the profile scale radius

# Keyword Arguments
- `dims...`: Optional dimensions for the output array (e.g., `100` for 100 samples, or `(10, 10)` for a 10×10 array)
- `rng`: Random number generator (default: `Random.default_rng()`)

# Returns
A vector or array of `SkyCoords` objects representing the positions of sampled points on the sky.

# Details
The function works by:
1. Sampling projected radii from the surface density profile of `d`
2. Converting physical radii to angular separations using the distance
3. Randomly distributing points uniformly in position angle (0 to 2π)
4. Using `SkyCoords.offset` to compute sky coordinates relative to the galaxy center

# Example
```julia
using GalaxyProfiles, SkyCoords
import Random

# Create a Plummer profile with scale radius 500 pc  
plummer = Plummer(1000.0, 0.5)

# Galaxy center at RA=270°, Dec=66° (converted to radians)
center = ICRSCoords(deg2rad(270.0), deg2rad(66.0))

# Distance to galaxy: 100 kpc
distance = 100.0

# Sample 1000 positions on the sky
coords = sample2D(plummer, center, distance, 1000)

# These coordinates are of the same type as the input `center`,
# in radians. We can separate out into RA and Dec and convert to deg2rad
# if needed.
ra, dec = rad2deg.(SkyCoords.lon.(coords)), rad2deg.(SkyCoords.lat.(coords))
```
"""
function sample2D(rng::AbstractRNG, d::AbstractMassProfile, center::AbstractSkyCoords,
                  distance)
    # Sample a single projected radius
    r_proj = sample2D_r(rng, d)
    # Sample random position angle uniformly from 0 to 2π
    pa = 2π * rand(rng)
    # Convert physical radius to angular separation (in radians)
    sep = atan(r_proj / distance)
    # Offset from center
    return offset(center, sep, pa)
end
sample2D(d::AbstractMassProfile, center::AbstractSkyCoords, distance) = 
    sample2D(default_rng(), d, center, distance)
function sample2D(rng::AbstractRNG, d::AbstractMassProfile, center::AbstractSkyCoords, 
                  distance, dims::Integer...)
    [sample2D(rng, d, center, distance) for _ in CartesianIndices(dims)]
end
function sample2D(d::AbstractMassProfile, center::AbstractSkyCoords, 
                  distance, dims::Integer...)
    sample2D(default_rng(), d, center, distance, dims...)
end
function sample2D(rng::AbstractRNG, d::AbstractMassProfile, center::AbstractSkyCoords,
                  distance::Real, dims::Dims)
    sample2D(rng, d, center, distance, dims...)
end
sample2D(d::AbstractMassProfile, center::AbstractSkyCoords, distance::Real, dims::Dims) =
    sample2D(default_rng(), d, center, distance, dims...)
    
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
