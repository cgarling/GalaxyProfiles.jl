module GalaxyProfiles

# definitions of functions to be overloaded in submodules
include("common.jl")
# includes the SurfaceDensities module
include("surface_densities/surface_densities.jl")
using .SurfaceDensities

export ExponentialDisk
export Σ

end # module
