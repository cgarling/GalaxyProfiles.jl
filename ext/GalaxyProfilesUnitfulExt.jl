module GalaxyProfilesUnitfulExt # Same name as file
# See https://docs.julialang.org/en/v1/manual/modules/#Submodules-and-relative-paths
# for relative module path conventions

if isdefined(Base, :get_extension)
    import GalaxyProfiles: AbstractMassProfile, AbstractDensity, ExponentialDisk, ExponentialDiskDHI, GeneralIsothermal, SIS, NFW, scale_radius, ρ, invρ, ∇ρ, ρmean, invρmean, Σ, ∇Σ, Σmean, invΣ, M, ∇M, invM, Mtot, Mproj, ∇Mproj, invMproj, dynamical_time, Vcirc, Vesc, Vmax, σr, σlos, Φ, ∇Φ, ∇∇Φ
    import Unitful as u
    import UnitfulAstro as ua
else
    # Up one module = ..
    import ..GalaxyProfiles: AbstractMassProfile, AbstractDensity, ExponentialDisk, ExponentialDiskDHI, GeneralIsothermal, SIS, NFW, scale_radius, ρ, invρ, ∇ρ, ρmean, invρmean, Σ, ∇Σ, Σmean, invΣ, M, ∇M, invM, Mtot, Mproj, ∇Mproj, invMproj, dynamical_time, Vcirc, Vesc, Vmax, σr, σlos, Φ, ∇Φ, ∇∇Φ
    import ..Unitful as u
    import ..UnitfulAstro as ua
end

# Define dimensionalities for dispatch
# 1*ua.Msun/ua.pc^2 isa GalaxyProfiles.SurfaceDensity
# ua.Msun/ua.pc^2 isa GalaxyProfiles.SurfaceDensityUnits
# Things like SurfaceDensityUnits are created automatically by @derived_dimension SurfaceDensity u.𝐌/u.𝐋^2
u.@derived_dimension SurfaceDensity u.𝐌/u.𝐋^2
u.@derived_dimension ∇ρdimension u.𝐌/u.𝐋^4
u.@derived_dimension Φdimension u.𝐋^2/u.𝐓^2
u.@derived_dimension ∇∇Φdimension u.𝐓^-2
u.@derived_dimension ∇Mdimension u.𝐌/u.𝐋
# u.@derived_dimension ∇Φdimension u.𝐋/u.𝐓^2  # this is just u.AccelerationUnits

# module to hold default units
module defaultunits

# Up two modules = ...
# import ...UnitfulAstro as ua
# import ...Unitful as u
# Actually, just import u and ua from one module above
import ..u # Not sure why this doesn't work with package extension
import ..ua

const time = u.yr # for timescales
const mass = ua.Msun
const ∇mass = ua.Msun / ua.kpc
const density = ua.Msun / ua.kpc^3
const ∇density = ua.Msun / ua.kpc^4
const surfacedensity = ua.Msun / ua.kpc^2
const length = ua.kpc
const velocity = u.km / u.s
const Φunit = u.km^2 / u.s^2
const ∇Φunit = u.km/u.s^2 # u.km^2 / u.s^2 / ua.kpc
const ∇∇Φunit = u.km/u.s^2/ua.kpc # u.km^2 / u.s^2 / ua.kpc^2

end # defaultunits module


# General functions for converting units
# Should check out how to have user-controllable default units, rather than the hard-coded values here
# convert_density(ρ0::u.Density) = u.ustrip(ua.Msun/ua.kpc^3, ρ0)
# convert_mass(M::u.Mass) = u.ustrip(ua.Msun, M)
# convert_surface_density(Σ0::SurfaceDensity) = u.ustrip(ua.Msun/ua.kpc^2, Σ0)
# convert_length(r::u.Length) = u.ustrip(ua.kpc, r)
homogenize_units(ρ::u.Density) = u.ustrip(defaultunits.density, ρ)
homogenize_units(M::u.Mass) = u.ustrip(defaultunits.mass, M)
homogenize_units(Σ::SurfaceDensity) = u.ustrip(defaultunits.surfacedensity, Σ)
homogenize_units(r::u.Length) = u.ustrip(defaultunits.length, r)
# Fallback for quantities with no default units defined
homogenize_units(x::u.Quantity) = throw(ArgumentError("No `homogenize_units` rule for input."))
# Fallback for non-unit numbers
homogenize_units(x::Real) = x

# Create Unitful constructors for our various composite types
ExponentialDisk(Σ0::SurfaceDensity, rs::u.Length) = ExponentialDisk(homogenize_units(Σ0), homogenize_units(rs))
function ExponentialDisk(rs::u.Length; M=nothing, Σ0=nothing)
    if isnothing(M)
        @assert !isnothing(Σ0)
        @assert Σ0 isa SurfaceDensity
        ExponentialDisk(Σ0, rs)
    else
        @assert M isa u.Mass
        # ExponentialDisk(M/(2π*rs^2),rs)
        ExponentialDisk(M / rs^2 / 2 / π, rs)
    end
end
ExponentialDiskDHI(DHI::u.Length, MHI::u.Mass, ΣDHI::SurfaceDensity=1*ua.Msun/ua.pc^2) = ExponentialDiskDHI(homogenize_units(DHI), homogenize_units(MHI), homogenize_units(ΣDHI))

GeneralIsothermal(ρ0::u.Density, rs::u.Length, α::Real) = GeneralIsothermal(homogenize_units(ρ0), homogenize_units(rs), α)
GeneralIsothermal(rs::u.Length, α::Real, M::u.Mass, Rmax::u.Length) = GeneralIsothermal(homogenize_units(rs), α, homogenize_units(M), homogenize_units(Rmax))

SIS(ρ0::u.Density, rs::u.Length) = SIS(homogenize_units(ρ0), homogenize_units(rs))
SIS(rs::u.Length, M::u.Mass, Rmax::u.Length) = SIS(homogenize_units(rs), homogenize_units(M), homogenize_units(Rmax))

NFW(ρ0::u.Density, rs::u.Length) = NFW(homogenize_units(ρ0), homogenize_units(rs))

#########################################################################################

scale_radius(uu::u.LengthUnits, d::AbstractMassProfile) = scale_radius(d) * defaultunits.length |> uu
ρ(d::AbstractDensity, r::u.Length) = ρ(d, homogenize_units(r)) * defaultunits.density
ρ(uu::u.DensityUnits, d::AbstractDensity, r::u.Length) = ρ(d, r) |> uu
ρ(uu::u.DensityUnits, d::AbstractDensity, r::Real) = ρ(d, r) * defaultunits.density |> uu
invρ(d::AbstractDensity, x::u.Density) = invρ(d, homogenize_units(x)) * defaultunits.length
invρ(uu::u.LengthUnits, d::AbstractDensity, x::u.Density) = invρ(d, x) |> uu
invρ(uu::u.LengthUnits, d::AbstractDensity, x::Real) = invρ(d, x) * defaultunits.length |> uu
∇ρ(d::AbstractDensity, r::u.Length) = ∇ρ(d, homogenize_units(r)) * defaultunits.∇density
∇ρ(uu::∇ρdimensionUnits, d::AbstractDensity, r::u.Length) = ∇ρ(d, r) |> uu
∇ρ(uu::∇ρdimensionUnits, d::AbstractDensity, r::Real) = ∇ρ(d, r) * defaultunits.∇density |> uu
ρmean(d::AbstractDensity, r::u.Length) = ρmean(d, homogenize_units(r)) * defaultunits.density
ρmean(uu::u.DensityUnits, d::AbstractDensity, r::u.Length) = ρmean(d, r) |> uu
ρmean(uu::u.DensityUnits, d::AbstractDensity, r::Real) = ρmean(d, r) * defaultunits.density |> uu
invρmean(d::AbstractDensity, x::u.Density) = invρmean(d, homogenize_units(x)) * defaultunits.length
invρmean(uu::u.LengthUnits, d::AbstractDensity, x::u.Density) = invρmean(d, x) |> uu
invρmean(uu::u.LengthUnits, d::AbstractDensity, x::Real) = invρmean(d, x) * defaultunits.length |> uu
Σ(d::AbstractMassProfile, r::u.Length) = Σ(d, homogenize_units(r)) * defaultunits.surfacedensity
Σ(uu::SurfaceDensityUnits, d::AbstractMassProfile, r::u.Length) = Σ(d, r) |> uu
Σ(uu::SurfaceDensityUnits, d::AbstractMassProfile, r::Real) = Σ(d, r) * defaultunits.surfacedensity |> uu
∇Σ(d::AbstractMassProfile, r::u.Length) = ∇Σ(d, homogenize_units(r)) * defaultunits.density
∇Σ(uu::u.DensityUnits, d::AbstractMassProfile, r::u.Length) = ∇Σ(d, r) |> uu
∇Σ(uu::u.DensityUnits, d::AbstractMassProfile, r::Real) = ∇Σ(d, r) * defaultunits.density |> uu
Σmean(d::AbstractMassProfile, r::u.Length) = Σmean(d, homogenize_units(r)) * defaultunits.surfacedensity
Σmean(uu::SurfaceDensityUnits, d::AbstractMassProfile, r::u.Length) = Σmean(d, r) |> uu
Σmean(uu::SurfaceDensityUnits, d::AbstractMassProfile, r::Real) = Σmean(d, r) * defaultunits.surfacedensity |> uu
invΣ(d::AbstractMassProfile, x::SurfaceDensity) = invΣ(d, homogenize_units(x)) * defaultunits.length
invΣ(uu::u.LengthUnits, d::AbstractMassProfile, x::SurfaceDensity) = invΣ(d, x) |> uu
invΣ(uu::u.LengthUnits, d::AbstractMassProfile, x::Real) = invΣ(d, x) * defaultunits.length |> uu
M(d::AbstractDensity, r::u.Length) = M(d, homogenize_units(r)) * defaultunits.mass
M(uu::u.MassUnits, d::AbstractDensity, r::u.Length) = M(d, r) |> uu
M(uu::u.MassUnits, d::AbstractDensity, r::Real) = M(d, r) * defaultunits.mass |> uu
∇M(d::AbstractDensity, r::u.Length) = ∇M(d, homogenize_units(r)) * defaultunits.∇mass
∇M(uu::∇MdimensionUnits, d::AbstractDensity, r::u.Length) = ∇M(d, r) |> uu
∇M(uu::∇MdimensionUnits, d::AbstractDensity, r::Real) = ∇M(d, r) * defaultunits.∇mass |> uu
invM(d::AbstractDensity, x::u.Mass) = invM(d, homogenize_units(x)) * defaultunits.length
invM(uu::u.LengthUnits, d::AbstractMassProfile, x::u.Mass) = invM(d, x) |> uu
invM(uu::u.LengthUnits, d::AbstractMassProfile, x::Real) = invM(d, x) * defaultunits.length |> uu
Mproj(d::AbstractDensity, r::u.Length) = Mproj(d, homogenize_units(r)) * defaultunits.mass
Mproj(uu::u.MassUnits, d::AbstractDensity, r::u.Length) = Mproj(d, r) |> uu
Mproj(uu::u.MassUnits, d::AbstractDensity, r::Real) = Mproj(d, r) * defaultunits.mass |> uu
∇Mproj(d::AbstractDensity, r::u.Length) = ∇Mproj(d, homogenize_units(r)) * defaultunits.∇mass
∇Mproj(uu::∇MdimensionUnits, d::AbstractDensity, r::u.Length) = ∇Mproj(d, r) |> uu
∇Mproj(uu::∇MdimensionUnits, d::AbstractDensity, r::Real) = ∇Mproj(d, r) * defaultunits.∇mass |> uu
invMproj(d::AbstractDensity, x::u.Mass) = invMproj(d, homogenize_units(x)) * defaultunits.length
invMproj(uu::u.LengthUnits, d::AbstractMassProfile, x::u.Mass) = invMproj(d, x) |> uu
invMproj(uu::u.LengthUnits, d::AbstractMassProfile, x::Real) = invMproj(d, x) * defaultunits.length |> uu
dynamical_time(d::AbstractMassProfile, r::u.Length) = dynamical_time(d, homogenize_units(r)) * defaultunits.time
dynamical_time(uu::u.TimeUnits, d::AbstractMassProfile, r::u.Length) = dynamical_time(d, r) |> uu
dynamical_time(uu::u.TimeUnits, d::AbstractMassProfile, r::Real) = dynamical_time(d, r) * defaultunits.time |> uu
Vcirc(d::AbstractDensity, r::u.Length) = Vcirc(d, homogenize_units(r)) * defaultunits.velocity
Vcirc(uu::u.VelocityUnits, d::AbstractDensity, r::u.Length) = Vcirc(d, r) |> uu
Vcirc(uu::u.VelocityUnits, d::AbstractDensity, r::Real) = Vcirc(d, r) * defaultunits.velocity |> uu
Vesc(d::AbstractDensity, r::u.Length) = Vesc(d, homogenize_units(r)) * defaultunits.velocity
Vesc(uu::u.VelocityUnits, d::AbstractDensity, r::u.Length) = Vesc(d, r) |> uu
Vesc(uu::u.VelocityUnits, d::AbstractDensity, r::Real) = Vesc(d, r) * defaultunits.velocity |> uu
# Vmax(uu::u.VelocityUnits, d::AbstractDensity) = ... returns a tuple of (Vmax, r (Vcirc == Vmax))
σr(d::AbstractDensity, r::u.Length, β) = σr(d, homogenize_units(r), β) * defaultunits.velocity
σr(uu::u.VelocityUnits, d::AbstractDensity, r::u.Length, β) = σr(d, r, β) |> uu
σr(uu::u.VelocityUnits, d::AbstractDensity, r::Real, β) = σr(d, r, β) * defaultunits.velocity |> uu
σlos(d::AbstractDensity, r::u.Length, β) = σlos(d, homogenize_units(r), β) * defaultunits.velocity
σlos(uu::u.VelocityUnits, d::AbstractDensity, r::u.Length, β) = σlos(d, r, β) |> uu
σlos(uu::u.VelocityUnits, d::AbstractDensity, r::Real, β) = σlos(d, r, β) * defaultunits.velocity |> uu
Φ(d::AbstractDensity, r::u.Length) = Φ(d, homogenize_units(r)) * defaultunits.Φunit
Φ(uu::ΦdimensionUnits, d::AbstractDensity, r::u.Length) = Φ(d, r) |> uu
Φ(uu::ΦdimensionUnits, d::AbstractDensity, r::Real) = Φ(d, r) * defaultunits.Φunit |> uu
∇Φ(d::AbstractDensity, r::u.Length) = ∇Φ(d, homogenize_units(r)) * defaultunits.∇Φunit
∇Φ(uu::u.AccelerationUnits, d::AbstractDensity, r::u.Length) = ∇Φ(d, r) |> uu
∇Φ(uu::u.AccelerationUnits, d::AbstractDensity, r::Real) = ∇Φ(d, r) * defaultunits.∇Φunit |> uu
∇∇Φ(d::AbstractDensity, r::u.Length) = ∇∇Φ(d, homogenize_units(r)) * defaultunits.∇∇Φunit
∇∇Φ(uu::∇∇ΦdimensionUnits, d::AbstractDensity, r::u.Length) = ∇∇Φ(d, r) |> uu
∇∇Φ(uu::∇∇ΦdimensionUnits, d::AbstractDensity, r::Real) = ∇∇Φ(d, r) * defaultunits.∇∇Φunit |> uu

#########################################################################################
# for f in (:ExponentialDisk, :GeneralIsothermal, :Plummer, :NFW) # quantities requiring 3D densities
#     @eval if hasmethod(scale_radius,($f,))
#         scale_radius(uu::u.LengthUnits,$f,args...;kws...) = scale_radius($f,args...;kws...) * defaultunits.length |> uu
#     end
#     @eval if hasmethod(ρ,($f,Real)) # check if this method is defined for the current type
#         @eval ρ(uu::u.DensityUnits,$f,args...;kws...) = ρ($f,args...;kws...) * defaultunits.density |> uu
#         @eval ρ(uu::u.DensityUnits,$f,r::u.Length,args...;kws...) = ρ($f,homogenize_units(r),args...;kws...) * defaultunits.density |> uu
#         @eval ρ($f,r::u.Length,args...;kws...) = ρ($f,homogenize_units(r),args...;kws...) * defaultunits.density
#     end
#     @eval if hasmethod(invρ,($f,Real)) # check if this method is defined for the current type
#         @eval invρ(uu::u.LengthUnits,$f,args...;kws...) = invρ($f,args...;kws...) * defaultunits.length |> uu
#         @eval invρ(uu::u.LengthUnits,$f,x::u.Density,args...;kws...) = invρ($f,homogenize_units(x),args...;kws...) * defaultunits.length |> uu
#         @eval invρ($f,x::u.Density,args...;kws...) = invρ($f,homogenize_units(x),args...;kws...) * defaultunits.length
#     end
#     @eval if hasmethod(∇ρ,($f,Real)) # check if this method is defined for the current type
#         @eval ∇ρ(uu::∇ρdimensionUnits,$f,args...;kws...) = ∇ρ($f,args...;kws...) * defaultunits.∇density |> uu
#         @eval ∇ρ(uu::∇ρdimensionUnits,$f,r::u.Length,args...;kws...) = ∇ρ($f,homogenize_units(r),args...;kws...) * defaultunits.∇density |> uu
#         @eval ∇ρ($f,r::u.Length,args...;kws...) = ∇ρ($f,homogenize_units(r),args...;kws...) * defaultunits.∇density
#     end
#     @eval if hasmethod(ρmean,($f,Real))
#         @eval ρmean(uu::u.DensityUnits,$f,args...;kws...) = ρmean($f,args...;kws...) * defaultunits.density |> uu
#         @eval ρmean(uu::u.DensityUnits,$f,r::u.Length,args...;kws...) = ρmean($f,homogenize_units(r),args...;kws...) * defaultunits.density |> uu
#         @eval ρmean($f,r::u.Length,args...;kws...) = ρmean($f,homogenize_units(r),args...;kws...) * defaultunits.density
#     end
#     @eval if hasmethod(invρmean,($f,Real)) # check if this method is defined for the current type
#         @eval invρmean(uu::u.LengthUnits,$f,args...;kws...) = invρmean($f,args...;kws...) * defaultunits.length |> uu
#         @eval invρmean(uu::u.LengthUnits,$f,x::u.Density,args...;kws...) = invρmean($f,homogenize_units(x),args...;kws...) * defaultunits.length |> uu
#         @eval invρmean($f,x::u.Density,args...;kws...) = invρmean($f,homogenize_units(x),args...;kws...) * defaultunits.length
#     end
#     @eval if hasmethod(Σ,($f,Real)) # check if this method is defined for the current type
#         @eval Σ(uu::SurfaceDensityUnits,$f,args...;kws...) = Σ($f,args...;kws...) * defaultunits.surfacedensity |> uu
#         @eval Σ(uu::SurfaceDensityUnits,$f,r::u.Length,args...;kws...) = Σ($f,homogenize_units(r),args...;kws...) * defaultunits.surfacedensity |> uu
#         @eval Σ($f,r::u.Length,args...;kws...) = Σ($f,homogenize_units(r),args...;kws...) * defaultunits.surfacedensity
#     end
#     @eval if hasmethod(∇Σ,($f,Real)) # check if this method is defined for the current type
#         @eval ∇Σ(uu::u.DensityUnits,$f,args...;kws...) = ∇Σ($f,args...;kws...) * defaultunits.density |> uu
#         @eval ∇Σ(uu::u.DensityUnits,$f,r::u.Length,args...;kws...) = ∇Σ($f,homogenize_units(r),args...;kws...) * defaultunits.density |> uu
#         @eval ∇Σ($f,r::u.Length,args...;kws...) = ∇Σ($f,homogenize_units(r),args...;kws...) * defaultunits.density
#     end
#     @eval if hasmethod(Σmean,($f,Real)) # check if this method is defined for the current type
#         @eval Σmean(uu::SurfaceDensityUnits,$f,args...;kws...) = Σmean($f,args...;kws...) * defaultunits.surfacedensity |> uu
#         @eval Σmean(uu::SurfaceDensityUnits,$f,r::u.Length,args...;kws...) = Σmean($f,homogenize_units(r),args...;kws...) * defaultunits.surfacedensity |> uu
#         @eval Σmean($f,r::u.Length,args...;kws...) = Σmean($f,homogenize_units(r),args...;kws...) * defaultunits.surfacedensity
#     end
#     @eval if hasmethod(invΣ,($f,Real)) # check if this method is defined for the current type
#         @eval invΣ(uu::u.LengthUnits,$f,args...;kws...) = invΣ($f,args...;kws...) * defaultunits.length |> uu
#         @eval invΣ(uu::u.LengthUnits,$f,x::SurfaceDensity,args...;kws...) = invΣ($f,homogenize_units(x),args...;kws...) * defaultunits.length |> uu
#         @eval invΣ($f,x::SurfaceDensity,args...;kws...) = invΣ($f,homogenize_units(x),args...;kws...) * defaultunits.length
#     end
#     @eval if hasmethod(M,($f,Real)) # check if this method is defined for the current type
#         @eval M(uu::u.MassUnits,$f,args...;kws...) = M($f,args...;kws...) * defaultunits.mass |> uu
#         @eval M(uu::u.MassUnits,$f,r::u.Length,args...;kws...) = M($f,homogenize_units(r),args...;kws...) * defaultunits.mass |> uu
#         @eval M($f,r::u.Length,args...;kws...) = M($f,homogenize_units(r),args...;kws...) * defaultunits.mass
#     end
#     @eval if hasmethod(∇M,($f,Real)) # check if this method is defined for the current type
#         @eval ∇M(uu::∇mdimensionUnits,$f,args...;kws...) = ∇M($f,args...;kws...) * defaultunits.∇mass |> uu
#         @eval ∇M(uu::∇mdimensionUnits,$f,r::u.Length,args...;kws...) = ∇M($f,homogenize_units(r),args...;kws...) * defaultunits.∇mass |> uu
#         @eval ∇M($f,r::u.Length,args...;kws...) = ∇M($f,homogenize_units(r),args...;kws...) * defaultunits.∇mass
#     end
#     @eval if hasmethod(invM,($f,Real)) # check if this method is defined for the current type
#         @eval invM(uu::u.LengthUnits,$f,args...;kws...) = invM($f,args...;kws...) * defaultunits.length |> uu
#         @eval invM(uu::u.LengthUnits,$f,x::u.Mass,args...;kws...) = invM($f,homogenize_units(x),args...;kws...) * defaultunits.length |> uu
#         @eval invM($f,x::u.Mass,args...;kws...) = invM($f,homogenize_units(x),args...;kws...) * defaultunits.length
#     end
#     @eval if hasmethod(Mtot,($f,)) # check if this method is defined for the current type
#         @eval Mtot(uu::u.MassUnits,$f,args...;kws...) = Mtot($f,args...;kws...) * defaultunits.mass |> uu
#     end
#     @eval if hasmethod(Mproj,($f,Real)) # check if this method is defined for the current type
#         @eval Mproj(uu::u.MassUnits,$f,args...;kws...) = Mproj($f,args...;kws...) * defaultunits.mass |> uu
#         @eval Mproj(uu::u.MassUnits,$f,r::u.Length,args...;kws...) = Mproj($f,homogenize_units(r),args...;kws...) * defaultunits.mass |> uu
#         @eval Mproj($f,r::u.Length,args...;kws...) = Mproj($f,homogenize_units(r),args...;kws...) * defaultunits.mass
#     end
#     @eval if hasmethod(∇Mproj,($f,Real)) # check if this method is defined for the current type
#         @eval ∇Mproj(uu::∇mdimensionUnits,$f,args...;kws...) = ∇Mproj($f,args...;kws...) * defaultunits.∇mass |> uu
#         @eval ∇Mproj(uu::∇mdimensionUnits,$f,r::u.Length,args...;kws...) = ∇Mproj($f,homogenize_units(r),args...;kws...) * defaultunits.∇mass |> uu
#         @eval ∇Mproj($f,r::u.Length,args...;kws...) = ∇Mproj($f,homogenize_units(r),args...;kws...) * defaultunits.∇mass
#     end
#     @eval if hasmethod(invMproj,($f,Real)) # check if this method is defined for the current type
#         @eval invMproj(uu::u.LengthUnits,$f,args...;kws...) = invMproj($f,args...;kws...) * defaultunits.length |> uu
#         @eval invMproj(uu::u.LengthUnits,$f,x::u.Mass,args...;kws...) = invMproj($f,homogenize_units(x),args...;kws...) * defaultunits.length |> uu
#         @eval invMproj($f,x::u.Mass,args...;kws...) = invMproj($f,homogenize_units(x),args...;kws...) * defaultunits.length
#     end
#     @eval if hasmethod(dynamical_time,($f,Real)) # check if this method is defined for the current type
#         @eval dynamical_time(uu::u.TimeUnits,$f,args...;kws...) = dynamical_time($f,args...;kws...) * defaultunits.time |> uu
#         @eval dynamical_time(uu::u.TimeUnits,$f,x::u.Length,args...;kws...) = dynamical_time($f,homogenize_units(x),args...;kws...) * defaultunits.time |> uu
#         @eval dynamical_time($f,x::u.Length,args...;kws...) = dynamical_time($f,homogenize_units(x),args...;kws...) * defaultunits.time
#     end
#     @eval if hasmethod(Vcirc,($f,Real))
#         @eval Vcirc(uu::u.VelocityUnits,$f,args...;kws...) = Vcirc($f,args...;kws...) * defaultunits.velocity |> uu
#         @eval Vcirc(uu::u.VelocityUnits,$f,r::u.Length,args...;kws...) = Vcirc($f,homogenize_units(r),args...;kws...) * defaultunits.velocity |> uu
#         @eval Vcirc($f,r::u.Length,args...;kws...) = Vcirc($f,homogenize_units(r),args...;kws...) * defaultunits.velocity       
#     end
#     @eval if hasmethod(Vesc,($f,Real))
#         @eval Vesc(uu::u.VelocityUnits,$f,args...;kws...) = Vesc($f,args...;kws...) * defaultunits.velocity |> uu
#         @eval Vesc(uu::u.VelocityUnits,$f,r::u.Length,args...;kws...) = Vesc($f,homogenize_units(r),args...;kws...) * defaultunits.velocity |> uu
#         @eval Vesc($f,r::u.Length,args...;kws...) = Vesc($f,homogenize_units(r),args...;kws...) * defaultunits.velocity       
#     end
#     @eval if hasmethod(Vmax,($f,Real))
#         @eval Vmax(uu::u.VelocityUnits,$f,args...;kws...) = Vmax($f,args...;kws...) * defaultunits.velocity |> uu
#         @eval Vmax(uu::u.VelocityUnits,$f,r::u.Length,args...;kws...) = Vmax($f,homogenize_units(r),args...;kws...) * defaultunits.velocity |> uu
#         @eval Vmax($f,r::u.Length,args...;kws...) = Vmax($f,homogenize_units(r),args...;kws...) * defaultunits.velocity
#     end
#     @eval if hasmethod(σr,($f,Real,Real))
#         @eval σr(uu::u.VelocityUnits,$f,args...;kws...) = σr($f,args...;kws...) * defaultunits.velocity |> uu
#         @eval σr(uu::u.VelocityUnits,$f,r::u.Length,args...;kws...) = σr($f,homogenize_units(r),args...;kws...) * defaultunits.velocity |> uu
#         @eval σr($f,r::u.Length,args...;kws...) = σr($f,homogenize_units(r),args...;kws...) * defaultunits.velocity
#     end
#     @eval if hasmethod(σlos,($f,Real,Real))
#         @eval σlos(uu::u.VelocityUnits,$f,args...;kws...) = σlos($f,args...;kws...) * defaultunits.velocity |> uu
#         @eval σlos(uu::u.VelocityUnits,$f,r::u.Length,args...;kws...) = σlos($f,homogenize_units(r),args...;kws...) * defaultunits.velocity |> uu
#         @eval σlos($f,r::u.Length,args...;kws...) = σlos($f,homogenize_units(r),args...;kws...) * defaultunits.velocity
#     end
#     @eval if hasmethod(Φ,($f,Real)) # check if this method is defined for the current type
#         @eval Φ(uu::ΦdimensionUnits,$f,args...;kws...) = Φ($f,args...;kws...) * defaultunits.Φunit |> uu
#         @eval Φ(uu::ΦdimensionUnits,$f,r::u.Length,args...;kws...) = Φ($f,homogenize_units(r),args...;kws...) * defaultunits.Φunit |> uu
#         @eval Φ($f,r::u.Length,args...;kws...) = Φ($f,homogenize_units(r),args...;kws...) * defaultunits.Φunit
#     end
#     @eval if hasmethod(∇Φ,($f,Real)) # check if this method is defined for the current type
#         @eval ∇Φ(uu::u.AccelerationUnits,$f,args...;kws...) = ∇Φ($f,args...;kws...) * defaultunits.∇Φunit |> uu
#         @eval ∇Φ(uu::u.AccelerationUnits,$f,r::u.Length,args...;kws...) = ∇Φ($f,homogenize_units(r),args...;kws...) * defaultunits.∇Φunit |> uu
#         @eval ∇Φ($f,r::u.Length,args...;kws...) = ∇Φ($f,homogenize_units(r),args...;kws...) * defaultunits.∇Φunit
#     end
#     @eval if hasmethod(∇∇Φ,($f,Real)) # check if this method is defined for the current type
#         @eval ∇∇Φ(uu::∇∇ΦdimensionUnits,$f,args...;kws...) = ∇∇Φ($f,args...;kws...) * defaultunits.∇∇Φunit |> uu
#         @eval ∇∇Φ(uu::∇∇ΦdimensionUnits,$f,r::u.Length,args...;kws...) = ∇∇Φ($f,homogenize_units(r),args...;kws...) * defaultunits.∇∇Φunit |> uu
#         @eval ∇∇Φ($f,r::u.Length,args...;kws...) = ∇∇Φ($f,homogenize_units(r),args...;kws...) * defaultunits.∇∇Φunit
#     end
# end

end # module
