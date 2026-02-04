module GalaxyProfilesSkyCoords

if isdefined(Base, :get_extension)
    import GalaxyProfiles: AbstractMassProfile, sample2D_r, sample2D_r!
    import SkyCoords
    import Random: AbstractRNG, default_rng, rand, rand!
else
    # For Julia < 1.9 without package extensions
    import ..GalaxyProfiles: AbstractMassProfile, sample2D_r, sample2D_r!
    import ..SkyCoords
    import Random: AbstractRNG, default_rng, rand, rand!
end

"""
    sample2D_skycoords(d::AbstractMassProfile, center::SkyCoords.AbstractSkyCoords, distance; 
                       dims..., rng=Random.default_rng())

Sample 2D projected radii from the mass profile `d` and convert them to sky coordinates.

# Arguments
- `d::AbstractMassProfile`: The mass profile to sample from
- `center::SkyCoords.AbstractSkyCoords`: The center position of the galaxy (e.g., ICRSCoords)
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

# Create a Plummer profile with scale radius 10 kpc  
plummer = Plummer(1000.0, 10.0)

# Galaxy center at RA=180°, Dec=30° (converted to radians)
center = ICRSCoords(deg2rad(180.0), deg2rad(30.0))

# Distance to galaxy: 10 Mpc = 10000 kpc
distance = 10000.0

# Sample 100 positions on the sky
coords = sample2D_skycoords(plummer, center, distance, 100)
```
"""
function sample2D_skycoords end

# Implementation for non-mutating forms with various dim arguments
function sample2D_skycoords(rng::AbstractRNG, d::AbstractMassProfile, center::SkyCoords.AbstractSkyCoords, 
                            distance::Real, dims::Vararg{Int, N}) where N
    # Sample projected radii in physical units
    r_proj = sample2D_r(rng, d, dims...)
    
    # Sample random position angles uniformly from 0 to 2π
    pas = 2π .* rand(rng, dims...)
    
    # Convert physical radii to angular separations (in radians)
    seps = r_proj ./ distance
    
    # Apply offsets for each sample
    return map((s, p) -> SkyCoords.offset(center, s, p), seps, pas)
end

# Convenience wrappers without explicit rng
function sample2D_skycoords(d::AbstractMassProfile, center::SkyCoords.AbstractSkyCoords, 
                            distance::Real, dims::Vararg{Int, N}) where N
    sample2D_skycoords(default_rng(), d, center, distance, dims...)
end

# Handle Dims{N} (tuple of ints) argument
function sample2D_skycoords(rng::AbstractRNG, d::AbstractMassProfile, center::SkyCoords.AbstractSkyCoords,
                            distance::Real, dims::Dims)
    sample2D_skycoords(rng, d, center, distance, dims...)
end

function sample2D_skycoords(d::AbstractMassProfile, center::SkyCoords.AbstractSkyCoords,
                            distance::Real, dims::Dims)
    sample2D_skycoords(default_rng(), d, center, distance, dims)
end

# Single sample versions
function sample2D_skycoords(rng::AbstractRNG, d::AbstractMassProfile, center::SkyCoords.AbstractSkyCoords,
                            distance::Real)
    # Sample a single projected radius
    r_proj = sample2D_r(rng, d)
    # Sample random position angle uniformly from 0 to 2π
    pa = 2π * rand(rng)
    # Convert physical radius to angular separation (in radians)
    sep = r_proj / distance
    # Offset from center
    return SkyCoords.offset(center, sep, pa)
end

function sample2D_skycoords(d::AbstractMassProfile, center::SkyCoords.AbstractSkyCoords,
                            distance::Real)
    sample2D_skycoords(default_rng(), d, center, distance)
end

# Export the function
export sample2D_skycoords

end # module
