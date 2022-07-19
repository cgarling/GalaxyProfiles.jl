module SurfaceDensities

import ..GalaxyProfiles:Σ

""" Abstract type for surface density profiles for which 3D quantities are not defined. Radii for instances of `AbstractSurfaceDensity` are always 2D. """
abstract type AbstractSurfaceDensity end

include("exponential_disk.jl")
# export types
export ExponentialDisk
# export methods
export Σ

end
