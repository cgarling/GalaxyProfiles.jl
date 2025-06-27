module GalaxyProfiles

using FastPow: @fastpow
using HypergeometricFunctions: _₂F₁ # For Plummer σr
using Random: AbstractRNG,rand,rand!,default_rng
using Roots: find_zero
using LambertW: lambertw
using SpecialFunctions: gamma, gamma_inc, gamma_inc_inv
using QuadGK: quadgk
# This symbol is only defined on Julia versions that support extensions
if !isdefined(Base, :get_extension)
    using Requires: @require
end

""" `AbstractMassProfile{T <: Real}`: abstract supertype for all mass profiles. """
abstract type AbstractMassProfile{T <: Real} end
Base.Broadcast.broadcastable(m::AbstractMassProfile) = Ref(m)

""" `AbstractDensity{T} <: AbstractMassProfile{T}`: abstract supertype for all 3D density profiles. """
abstract type AbstractDensity{T} <: AbstractMassProfile{T} end

""" Abstract type `(<:AbstractMassProfile)` for surface density profiles for which 3D quantities are not defined. Radii for instances of `AbstractSurfaceDensity` are always 2D. """
abstract type AbstractSurfaceDensity{T} <: AbstractMassProfile{T} end

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
    @static if !isdefined(Base, :get_extension)
        @require Unitful="1986cc42-f94f-5a68-af5c-568840ba703d" begin
            @require UnitfulAstro="6112ee07-acf9-5e0f-b108-d242c714bf9f" include("../ext/GalaxyProfilesUnitfulExt.jl")
        end
    end
end

export ExponentialDisk, ExponentialDiskDHI, Sersic  # Exports from surface_densities/*
export GeneralIsothermal, SIS, NFW, CoreNFW, CoreNFWGalaxy, Plummer # Exports from densities/*
export params, scale_radius, ρ, ρmean, invρmean, ∇ρ, invρ, Σ, invΣ, ∇Σ, Σmean, M, ∇M, invM, Mtot, Mproj, ∇Mproj, invMproj, dynamical_time, cdf2D, cdf3D, ccdf2D, ccdf3D, quantile2D, quantile3D, cquantile2D, cquantile3D, Vcirc, Vesc, Vmax, σr, σlos, Φ, ∇Φ, ∇∇Φ  # Exports from common.jl
export sample2D_r!, sample3D_r!, sample2D_r, sample3D_r # Exports from generic_rand.jl.
end # module
