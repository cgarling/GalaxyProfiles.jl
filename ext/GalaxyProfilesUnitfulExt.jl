module GalaxyProfilesUnitfulExt # Same name as file
# See https://docs.julialang.org/en/v1/manual/modules/#Submodules-and-relative-paths
# for relative module path conventions

if isdefined(Base, :get_extension)
    import GalaxyProfiles: AbstractMassProfile, AbstractDensity, ExponentialDisk, ExponentialDiskDHI, Sersic, GeneralIsothermal, SIS, NFW, CoreNFW, Plummer, scale_radius, ρ, invρ, ∇ρ, ρmean, invρmean, Σ, ∇Σ, Σmean, invΣ, M, ∇M, invM, Mtot, Mproj, ∇Mproj, invMproj, dynamical_time, Vcirc, Vesc, Vmax, σr, σlos, Φ, ∇Φ, ∇∇Φ
    import Unitful as u
    import UnitfulAstro as ua
else # For Julia < 1.9 without package extensions
    # Up one module = ..
    import ..GalaxyProfiles: AbstractMassProfile, AbstractDensity, ExponentialDisk, ExponentialDiskDHI, Sersic, GeneralIsothermal, SIS, NFW, CoreNFW, Plummer, scale_radius, ρ, invρ, ∇ρ, ρmean, invρmean, Σ, ∇Σ, Σmean, invΣ, M, ∇M, invM, Mtot, Mproj, ∇Mproj, invMproj, dynamical_time, Vcirc, Vesc, Vmax, σr, σlos, Φ, ∇Φ, ∇∇Φ
    import ..Unitful as u
    import ..UnitfulAstro as ua
end

# Define dimensions for dispatch
# 1*ua.Msun/ua.pc^2 isa GalaxyProfiles.SurfaceDensity
# ua.Msun/ua.pc^2 isa GalaxyProfiles.SurfaceDensityUnits
# Things like SurfaceDensityUnits are created automatically by @derived_dimension SurfaceDensity u.𝐌/u.𝐋^2
u.@derived_dimension SurfaceDensity u.𝐌/u.𝐋^2
u.@derived_dimension ∇ρdimension u.𝐌/u.𝐋^4
u.@derived_dimension Φdimension u.𝐋^2/u.𝐓^2
u.@derived_dimension ∇∇Φdimension u.𝐓^-2
u.@derived_dimension ∇Mdimension u.𝐌/u.𝐋
# u.@derived_dimension ∇Φdimension u.𝐋/u.𝐓^2  # this is just u.AccelerationUnits

#########################################################################################

module defaultunits # module to hold default units

# Up two modules = ...
# import ...UnitfulAstro as ua
# import ...Unitful as u
# Actually, just import u and ua from one module above
import ..u
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

#########################################################################################

# General functions for converting units
homogenize_units(ρ::u.Density) = u.ustrip(defaultunits.density, ρ)
homogenize_units(M::u.Mass) = u.ustrip(defaultunits.mass, M)
homogenize_units(Σ::SurfaceDensity) = u.ustrip(defaultunits.surfacedensity, Σ)
homogenize_units(r::u.Length) = u.ustrip(defaultunits.length, r)
# Fallback for quantities with no default units defined
homogenize_units(x::u.Quantity) = throw(ArgumentError("No `homogenize_units` rule for input."))
# Fallback for non-unit numbers
homogenize_units(x::Real) = x

#########################################################################################

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

Sersic(Σ0::SurfaceDensity, r_e::u.Length, n::Real, q::Real) = Sersic(homogenize_units(Σ0), homogenize_units(r_e), n, q)

GeneralIsothermal(ρ0::u.Density, rs::u.Length, α::Real) = GeneralIsothermal(homogenize_units(ρ0), homogenize_units(rs), α)
GeneralIsothermal(rs::u.Length, α::Real, M::u.Mass, Rmax::u.Length) = GeneralIsothermal(homogenize_units(rs), α, homogenize_units(M), homogenize_units(Rmax))

SIS(ρ0::u.Density, rs::u.Length) = SIS(homogenize_units(ρ0), homogenize_units(rs))
SIS(rs::u.Length, M::u.Mass, Rmax::u.Length) = SIS(homogenize_units(rs), homogenize_units(M), homogenize_units(Rmax))

NFW(ρ0::u.Density, rs::u.Length) = NFW(homogenize_units(ρ0), homogenize_units(rs))
# Plummer constructor removed - units are now preserved in the struct directly
#Plummer(M::u.Mass, a::u.Length) = Plummer(homogenize_units(M), homogenize_units(a))
CoreNFW(ρ0::u.Density, rs::u.Length, rc::u.Length, n::Real) =
    CoreNFW(homogenize_units(ρ0), homogenize_units(rs), homogenize_units(rc), n)

#########################################################################################

# Common methods supporting Unitful arguments
scale_radius(uu::u.LengthUnits, d::AbstractMassProfile) = scale_radius(d) * defaultunits.length |> uu
# ρ methods removed - units are now handled naturally in the core methods
# ρ(d::AbstractDensity, r::u.Length) = ρ(d, homogenize_units(r)) * defaultunits.density
# ρ(uu::u.DensityUnits, d::AbstractDensity, r::u.Length) = ρ(d, r) |> uu
# ρ(uu::u.DensityUnits, d::AbstractDensity, r::Real) = ρ(d, r) * defaultunits.density |> uu
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
Mtot(uu::u.MassUnits, d::AbstractMassProfile) = Mtot(d) * defaultunits.mass |> uu
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

end # module
