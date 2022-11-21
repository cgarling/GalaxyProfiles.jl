module GalaxyProfiles

import Random: AbstractRNG,rand,rand!,default_rng
import Roots: find_zero
import LambertW: lambertw
import SpecialFunctions: gamma
import QuadGK: quadgk
import Requires: @require

""" Supertype for all mass profiles. """
abstract type AbstractMassProfile end
Base.Broadcast.broadcastable(m::AbstractMassProfile) = Ref(m)

""" Abstract type `(<:AbstractMassProfile)` for 3D density profiles. """
abstract type AbstractDensity <: AbstractMassProfile end

""" Abstract type `(<:AbstractMassProfile)` for surface density profiles for which 3D quantities are not defined. Radii for instances of `AbstractSurfaceDensity` are always 2D. """
abstract type AbstractSurfaceDensity <: AbstractMassProfile end

# mutable struct UnitDefaults{T<:u.FreeUnits{S,u.𝐋,nothing} where S, V<:u.FreeUnits{Z,u.𝐌,nothing} where Z}
#     length::T
#     mass::V
# end
# const DefaultUnits = Dict("length"=>ua.kpc,"mass"=>ua.Msun)

# definitions of functions to be overloaded in submodules
include("common.jl")
include("surface_densities/surface_densities.jl")
include("densities/densities.jl")

# generic random sampling
include("generic_rand.jl")

function __init__()
    @require Unitful="1986cc42-f94f-5a68-af5c-568840ba703d" begin
        @require UnitfulAstro="6112ee07-acf9-5e0f-b108-d242c714bf9f" include("units.jl")
    end
end

export ExponentialDisk, ExponentialDiskDHI  # Exports from surface_densities/*
export GeneralIsothermal, SIS, NFW, Plummer # Exports from densities/*
export params, scale_radius, ρ, ρmean, invρmean, ∇ρ, invρ, Σ, invΣ, ∇Σ, Σmean, M, ∇M, invM, Mtot, Mproj, ∇Mproj, invMproj, cdf2D, cdf3D, ccdf2D, ccdf3D, quantile2D, quantile3D, cquantile2D, cquantile3D, Vcirc, Vesc, Vmax, Φ, ∇Φ, ∇∇Φ  # Exports from common.jl
export sample2D_r!, sample3D_r!, sample2D_r, sample3D_r # Exports from generic_rand.jl.
end # module
