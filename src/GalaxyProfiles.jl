module GalaxyProfiles

import Random:AbstractRNG,rand,rand!,default_rng
import Roots:find_zero
import LambertW:lambertw
import SpecialFunctions:gamma

""" Supertype for all mass profiles. """
abstract type AbstractMassProfile end
Base.Broadcast.broadcastable(m::AbstractMassProfile) = Ref(m)

""" Abstract type `(<:AbstractMassProfile)` for 3D density profiles. """
abstract type AbstractDensity <: AbstractMassProfile end

""" Abstract type `(<:AbstractMassProfile)` for surface density profiles for which 3D quantities are not defined. Radii for instances of `AbstractSurfaceDensity` are always 2D. """
abstract type AbstractSurfaceDensity <: AbstractMassProfile end

# definitions of functions to be overloaded in submodules
include("common.jl")
include("surface_densities/surface_densities.jl")
include("densities/densities.jl")

# generic random sampling
include("generic_rand.jl")

export ExponentialDisk
export GeneralIsothermal, SIS
export params, scale_radius, ρ, dρ_dr, invρ, Σ, invΣ, dΣ_dr, M, dM_dr, invM, Mtot, cdf, ccdf, quantile, cquantile, rand, rand!

end # module
