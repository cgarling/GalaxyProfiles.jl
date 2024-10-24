import .UnitfulAstro as ua
import .Unitful as u

# Define dimensionalities for dispatch
# 1*ua.Msun/ua.pc^2 isa GalaxyProfiles.SurfaceDensity
# ua.Msun/ua.pc^2 isa GalaxyProfiles.SurfaceDensityUnits
# This ...Units is created automatically by @derived_dimension.
u.@derived_dimension SurfaceDensity u.𝐌/u.𝐋^2
u.@derived_dimension ∇ρdimension u.𝐌/u.𝐋^4
u.@derived_dimension Φdimension u.𝐋^2/u.𝐓^2
u.@derived_dimension ∇∇Φdimension u.𝐓^-2
u.@derived_dimension ∇mdimension u.𝐌/u.𝐋
# u.@derived_dimension ∇Φdimension u.𝐋/u.𝐓^2  # this is just u.AccelerationUnits

# module to hold default units
module defaultunits
import ..UnitfulAstro as ua
import ..Unitful as u

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
end

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
        ExponentialDisk(M / rs^2 / 2 / π,rs)
    end
end
ExponentialDiskDHI(DHI::Unitful.Length, MHI::Unitful.Mass, ΣDHI::SurfaceDensity=1*ua.Msun/ua.pc^2) = ExponentialDiskDHI(homogenize_units(DHI), homogenize_units(MHI), homogenize_units(ΣDHI))

GeneralIsothermal(ρ0::u.Density, rs::u.Length, α::Real) = GeneralIsothermal(homogenize_units(ρ0), homogenize_units(rs), α)
GeneralIsothermal(rs::u.Length, α::Real, M::u.Mass, Rmax::u.Length) = GeneralIsothermal(homogenize_units(rs), α, homogenize_units(M), homogenize_units(Rmax))

SIS(ρ0::u.Density, rs::u.Length) = SIS(homogenize_units(ρ0), homogenize_units(rs))
SIS(rs::u.Length, M::u.Mass, Rmax::u.Length) = SIS(homogenize_units(rs), homogenize_units(M), homogenize_units(Rmax))

NFW(ρ0::u.Density, rs::u.Length) = NFW(homogenize_units(ρ0), homogenize_units(rs))

#########################################################################################
# Easily select a different unit
# This is pretty annoying but it's nice to have I guess ...
# This also has hard-coded default units, which I should probably specify somewhere explicitly
# Switch this to have loops over Types and over Functions

# for f in (:ExponentialDisk,:GeneralIsothermal) # common quantities for both 3D densities and 2D SurfaceDensities
#     @eval Σ(uu::u.Unitlike,$f,args...;kws...) = Σ($f,args...;kws...) * u.ustrip(1*defaultunits.surfacedensity |> uu)
#     @eval Σ(uu::u.Unitlike,$f,r::u.Length,args...;kws...) = Σ($f,u.ustrip(r|>defaultunits.length),args...;kws...) * u.ustrip(1*defaultunits.surfacedensity |> uu)
#     @eval Σ($f,r::u.Length,args...;kws...) = Σ($f,u.ustrip(r|>defaultunits.length),args...;kws...)
#     @eval ∇Σ(uu::u.Unitlike,$f,args...;kws...) = ∇Σ($f,args...;kws...) * u.ustrip(1*defaultunits.density |> uu)
#     @eval ∇Σ(uu::u.Unitlike,$f,r::u.Length,args...;kws...) = ∇Σ($f,u.ustrip(r|>defaultunits.length),args...;kws...) * u.ustrip(1*defaultunits.density |> uu)
#     @eval ∇Σ($f,r::u.Length,args...;kws...) = ∇Σ($f,u.ustrip(r|>defaultunits.length),args...;kws...)
# end
# for f in (:GeneralIsothermal,) # quantities requiring 3D densities
#     @eval Φ(uu::u.Unitlike,$f,args...;kws...) = Φ($f,args...;kws...) * u.ustrip(1*defaultunits.Φunit |> uu)
#     @eval Φ(uu::u.Unitlike,$f,r::u.Length,args...;kws...) = Φ($f,u.ustrip(r|>defaultunit.length),args...;kws...) * u.ustrip(1*defaultunits.Φunit |> uu)
#     @eval Φ($f,r::u.Length,args...;kws...) = Φ($f,u.ustrip(r|>defaultunit.length),args...;kws...)
#     @eval ∇Φ(uu::u.Unitlike,$f,args...;kws...) = ∇Φ($f,args...;kws...) * u.ustrip(1*defaultunits.∇Φunit |> uu)
#     @eval ∇Φ(uu::u.Unitlike,$f,r::u.Length,args...;kws...) = ∇Φ($f,u.ustrip(r|>defaultunits.length),args...;kws...) * u.ustrip(1*defaultunits.∇Φunit |> uu)
#     @eval ∇Φ($f,r::u.Length,args...;kws...) = ∇Φ($f,u.ustrip(r|>defaultunits.length),args...;kws...)
# end
# return quantities with units on them when requesting a different unit f(uu...), or providing a unit on an argument to the method f(d,[kpc])
# for f in (:ExponentialDisk,:GeneralIsothermal) # common quantities for both 3D densities and 2D SurfaceDensities
#     @eval Σ(uu::SurfaceDensityUnits,$f,args...;kws...) = Σ($f,args...;kws...) * defaultunits.surfacedensity |> uu
#     @eval Σ(uu::SurfaceDensityUnits,$f,r::u.Length,args...;kws...) = Σ($f,homogenize_units(r),args...;kws...) * defaultunits.surfacedensity |> uu
#     @eval Σ($f,r::u.Length,args...;kws...) = Σ($f,homogenize_units(r),args...;kws...) * defaultunits.surfacedensity
#     @eval ∇Σ(uu::u.DensityUnits,$f,args...;kws...) = ∇Σ($f,args...;kws...) * defaultunits.density |> uu
#     @eval ∇Σ(uu::u.DensityUnits,$f,r::u.Length,args...;kws...) = ∇Σ($f,homogenize_units(r),args...;kws...) * defaultunits.density |> uu
#     @eval ∇Σ($f,r::u.LengthUnits,args...;kws...) = ∇Σ($f,homogenize_units(r),args...;kws...) * defaultunits.density
#     # Σmean
#     @eval invΣ(uu::u.LengthUnits,$f,args...;kws...) = invΣ($f,args...;kws...) * defaultunits.length |> uu
#     @eval invΣ(uu::u.LengthUnits,$f,x::SurfaceDensity,args...;kws...) = invΣ($f,homogenize_units(x),args...;kws...) * defaultunits.length |> uu
#     @eval invΣ($f,x::SurfaceDensity,args...;kws...) = invΣ($f,homogenize_units(x),args...;kws...) * defaultunits.length
# end

# These definitions should allow MOST AbstractMassProfiles to have their unitful methods auto-generated.
# However, the args... are pointless in this context because we are specifying a specific signature with hasmethod.
# If the same methods are defined on other types which take more than one argument (or a single argument that is not a Real)
# then additional hasmethod branches will need to be defined.
for f in (:ExponentialDisk, :GeneralIsothermal, :Plummer, :NFW) # quantities requiring 3D densities
    @eval if hasmethod(scale_radius,($f,))
        scale_radius(uu::u.LengthUnits,$f,args...;kws...) = scale_radius($f,args...;kws...) * defaultunits.length |> uu
    end
    @eval if hasmethod(ρ,($f,Real)) # check if this method is defined for the current type
        @eval ρ(uu::u.DensityUnits,$f,args...;kws...) = ρ($f,args...;kws...) * defaultunits.density |> uu
        @eval ρ(uu::u.DensityUnits,$f,r::u.Length,args...;kws...) = ρ($f,homogenize_units(r),args...;kws...) * defaultunits.density |> uu
        @eval ρ($f,r::u.Length,args...;kws...) = ρ($f,homogenize_units(r),args...;kws...) * defaultunits.density
    end
    @eval if hasmethod(invρ,($f,Real)) # check if this method is defined for the current type
        @eval invρ(uu::u.LengthUnits,$f,args...;kws...) = invρ($f,args...;kws...) * defaultunits.length |> uu
        @eval invρ(uu::u.LengthUnits,$f,x::u.Density,args...;kws...) = invρ($f,homogenize_units(x),args...;kws...) * defaultunits.length |> uu
        @eval invρ($f,x::u.Density,args...;kws...) = invρ($f,homogenize_units(x),args...;kws...) * defaultunits.length
    end
    @eval if hasmethod(∇ρ,($f,Real)) # check if this method is defined for the current type
        @eval ∇ρ(uu::∇ρdimensionUnits,$f,args...;kws...) = ∇ρ($f,args...;kws...) * defaultunits.∇density |> uu
        @eval ∇ρ(uu::∇ρdimensionUnits,$f,r::u.Length,args...;kws...) = ∇ρ($f,homogenize_units(r),args...;kws...) * defaultunits.∇density |> uu
        @eval ∇ρ($f,r::u.Length,args...;kws...) = ∇ρ($f,homogenize_units(r),args...;kws...) * defaultunits.∇density
    end
    @eval if hasmethod(ρmean,($f,Real))
        @eval ρmean(uu::u.DensityUnits,$f,args...;kws...) = ρmean($f,args...;kws...) * defaultunits.density |> uu
        @eval ρmean(uu::u.DensityUnits,$f,r::u.Length,args...;kws...) = ρmean($f,homogenize_units(r),args...;kws...) * defaultunits.density |> uu
        @eval ρmean($f,r::u.Length,args...;kws...) = ρmean($f,homogenize_units(r),args...;kws...) * defaultunits.density
    end
    @eval if hasmethod(invρmean,($f,Real)) # check if this method is defined for the current type
        @eval invρmean(uu::u.LengthUnits,$f,args...;kws...) = invρmean($f,args...;kws...) * defaultunits.length |> uu
        @eval invρmean(uu::u.LengthUnits,$f,x::u.Density,args...;kws...) = invρmean($f,homogenize_units(x),args...;kws...) * defaultunits.length |> uu
        @eval invρmean($f,x::u.Density,args...;kws...) = invρmean($f,homogenize_units(x),args...;kws...) * defaultunits.length
    end
    @eval if hasmethod(Σ,($f,Real)) # check if this method is defined for the current type
        @eval Σ(uu::SurfaceDensityUnits,$f,args...;kws...) = Σ($f,args...;kws...) * defaultunits.surfacedensity |> uu
        @eval Σ(uu::SurfaceDensityUnits,$f,r::u.Length,args...;kws...) = Σ($f,homogenize_units(r),args...;kws...) * defaultunits.surfacedensity |> uu
        @eval Σ($f,r::u.Length,args...;kws...) = Σ($f,homogenize_units(r),args...;kws...) * defaultunits.surfacedensity
    end
    @eval if hasmethod(∇Σ,($f,Real)) # check if this method is defined for the current type
        @eval ∇Σ(uu::u.DensityUnits,$f,args...;kws...) = ∇Σ($f,args...;kws...) * defaultunits.density |> uu
        @eval ∇Σ(uu::u.DensityUnits,$f,r::u.Length,args...;kws...) = ∇Σ($f,homogenize_units(r),args...;kws...) * defaultunits.density |> uu
        @eval ∇Σ($f,r::u.Length,args...;kws...) = ∇Σ($f,homogenize_units(r),args...;kws...) * defaultunits.density
    end
    @eval if hasmethod(Σmean,($f,Real)) # check if this method is defined for the current type
        @eval Σmean(uu::SurfaceDensityUnits,$f,args...;kws...) = Σmean($f,args...;kws...) * defaultunits.surfacedensity |> uu
        @eval Σmean(uu::SurfaceDensityUnits,$f,r::u.Length,args...;kws...) = Σmean($f,homogenize_units(r),args...;kws...) * defaultunits.surfacedensity |> uu
        @eval Σmean($f,r::u.Length,args...;kws...) = Σmean($f,homogenize_units(r),args...;kws...) * defaultunits.surfacedensity
    end
    @eval if hasmethod(invΣ,($f,Real)) # check if this method is defined for the current type
        @eval invΣ(uu::u.LengthUnits,$f,args...;kws...) = invΣ($f,args...;kws...) * defaultunits.length |> uu
        @eval invΣ(uu::u.LengthUnits,$f,x::SurfaceDensity,args...;kws...) = invΣ($f,homogenize_units(x),args...;kws...) * defaultunits.length |> uu
        @eval invΣ($f,x::SurfaceDensity,args...;kws...) = invΣ($f,homogenize_units(x),args...;kws...) * defaultunits.length
    end
    @eval if hasmethod(M,($f,Real)) # check if this method is defined for the current type
        @eval M(uu::u.MassUnits,$f,args...;kws...) = M($f,args...;kws...) * defaultunits.mass |> uu
        @eval M(uu::u.MassUnits,$f,r::u.Length,args...;kws...) = M($f,homogenize_units(r),args...;kws...) * defaultunits.mass |> uu
        @eval M($f,r::u.Length,args...;kws...) = M($f,homogenize_units(r),args...;kws...) * defaultunits.mass
    end
    @eval if hasmethod(∇M,($f,Real)) # check if this method is defined for the current type
        @eval ∇M(uu::∇mdimensionUnits,$f,args...;kws...) = ∇M($f,args...;kws...) * defaultunits.∇mass |> uu
        @eval ∇M(uu::∇mdimensionUnits,$f,r::u.Length,args...;kws...) = ∇M($f,homogenize_units(r),args...;kws...) * defaultunits.∇mass |> uu
        @eval ∇M($f,r::u.Length,args...;kws...) = ∇M($f,homogenize_units(r),args...;kws...) * defaultunits.∇mass
    end
    @eval if hasmethod(invM,($f,Real)) # check if this method is defined for the current type
        @eval invM(uu::u.LengthUnits,$f,args...;kws...) = invM($f,args...;kws...) * defaultunits.length |> uu
        @eval invM(uu::u.LengthUnits,$f,x::u.Mass,args...;kws...) = invM($f,homogenize_units(x),args...;kws...) * defaultunits.length |> uu
        @eval invM($f,x::u.Mass,args...;kws...) = invM($f,homogenize_units(x),args...;kws...) * defaultunits.length
    end
    @eval if hasmethod(Mtot,($f,)) # check if this method is defined for the current type
        @eval Mtot(uu::u.MassUnits,$f,args...;kws...) = Mtot($f,args...;kws...) * defaultunits.mass |> uu
    end
    @eval if hasmethod(Mproj,($f,Real)) # check if this method is defined for the current type
        @eval Mproj(uu::u.MassUnits,$f,args...;kws...) = Mproj($f,args...;kws...) * defaultunits.mass |> uu
        @eval Mproj(uu::u.MassUnits,$f,r::u.Length,args...;kws...) = Mproj($f,homogenize_units(r),args...;kws...) * defaultunits.mass |> uu
        @eval Mproj($f,r::u.Length,args...;kws...) = Mproj($f,homogenize_units(r),args...;kws...) * defaultunits.mass
    end
    @eval if hasmethod(∇Mproj,($f,Real)) # check if this method is defined for the current type
        @eval ∇Mproj(uu::∇mdimensionUnits,$f,args...;kws...) = ∇Mproj($f,args...;kws...) * defaultunits.∇mass |> uu
        @eval ∇Mproj(uu::∇mdimensionUnits,$f,r::u.Length,args...;kws...) = ∇Mproj($f,homogenize_units(r),args...;kws...) * defaultunits.∇mass |> uu
        @eval ∇Mproj($f,r::u.Length,args...;kws...) = ∇Mproj($f,homogenize_units(r),args...;kws...) * defaultunits.∇mass
    end
    @eval if hasmethod(invMproj,($f,Real)) # check if this method is defined for the current type
        @eval invMproj(uu::u.LengthUnits,$f,args...;kws...) = invMproj($f,args...;kws...) * defaultunits.length |> uu
        @eval invMproj(uu::u.LengthUnits,$f,x::u.Mass,args...;kws...) = invMproj($f,homogenize_units(x),args...;kws...) * defaultunits.length |> uu
        @eval invMproj($f,x::u.Mass,args...;kws...) = invMproj($f,homogenize_units(x),args...;kws...) * defaultunits.length
    end
    @eval if hasmethod(dynamical_time,($f,Real)) # check if this method is defined for the current type
        @eval dynamical_time(uu::u.TimeUnits,$f,args...;kws...) = dynamical_time($f,args...;kws...) * defaultunits.time |> uu
        @eval dynamical_time(uu::u.TimeUnits,$f,x::u.Length,args...;kws...) = dynamical_time($f,homogenize_units(x),args...;kws...) * defaultunits.time |> uu
        @eval dynamical_time($f,x::u.Length,args...;kws...) = dynamical_time($f,homogenize_units(x),args...;kws...) * defaultunits.time
    end
    @eval if hasmethod(Vcirc,($f,Real))
        @eval Vcirc(uu::u.VelocityUnits,$f,args...;kws...) = Vcirc($f,args...;kws...) * defaultunits.velocity |> uu
        @eval Vcirc(uu::u.VelocityUnits,$f,r::u.Length,args...;kws...) = Vcirc($f,homogenize_units(r),args...;kws...) * defaultunits.velocity |> uu
        @eval Vcirc($f,r::u.Length,args...;kws...) = Vcirc($f,homogenize_units(r),args...;kws...) * defaultunits.velocity       
    end
    @eval if hasmethod(Vesc,($f,Real))
        @eval Vesc(uu::u.VelocityUnits,$f,args...;kws...) = Vesc($f,args...;kws...) * defaultunits.velocity |> uu
        @eval Vesc(uu::u.VelocityUnits,$f,r::u.Length,args...;kws...) = Vesc($f,homogenize_units(r),args...;kws...) * defaultunits.velocity |> uu
        @eval Vesc($f,r::u.Length,args...;kws...) = Vesc($f,homogenize_units(r),args...;kws...) * defaultunits.velocity       
    end
    @eval if hasmethod(Vmax,($f,Real))
        @eval Vmax(uu::u.VelocityUnits,$f,args...;kws...) = Vmax($f,args...;kws...) * defaultunits.velocity |> uu
        @eval Vmax(uu::u.VelocityUnits,$f,r::u.Length,args...;kws...) = Vmax($f,homogenize_units(r),args...;kws...) * defaultunits.velocity |> uu
        @eval Vmax($f,r::u.Length,args...;kws...) = Vmax($f,homogenize_units(r),args...;kws...) * defaultunits.velocity
    end
    @eval if hasmethod(σr,($f,Real,Real))
        @eval σr(uu::u.VelocityUnits,$f,args...;kws...) = σr($f,args...;kws...) * defaultunits.velocity |> uu
        @eval σr(uu::u.VelocityUnits,$f,r::u.Length,args...;kws...) = σr($f,homogenize_units(r),args...;kws...) * defaultunits.velocity |> uu
        @eval σr($f,r::u.Length,args...;kws...) = σr($f,homogenize_units(r),args...;kws...) * defaultunits.velocity
    end
    @eval if hasmethod(σlos,($f,Real,Real))
        @eval σlos(uu::u.VelocityUnits,$f,args...;kws...) = σlos($f,args...;kws...) * defaultunits.velocity |> uu
        @eval σlos(uu::u.VelocityUnits,$f,r::u.Length,args...;kws...) = σlos($f,homogenize_units(r),args...;kws...) * defaultunits.velocity |> uu
        @eval σlos($f,r::u.Length,args...;kws...) = σlos($f,homogenize_units(r),args...;kws...) * defaultunits.velocity
    end
    @eval if hasmethod(Φ,($f,Real)) # check if this method is defined for the current type
        @eval Φ(uu::ΦdimensionUnits,$f,args...;kws...) = Φ($f,args...;kws...) * defaultunits.Φunit |> uu
        @eval Φ(uu::ΦdimensionUnits,$f,r::u.Length,args...;kws...) = Φ($f,homogenize_units(r),args...;kws...) * defaultunits.Φunit |> uu
        @eval Φ($f,r::u.Length,args...;kws...) = Φ($f,homogenize_units(r),args...;kws...) * defaultunits.Φunit
    end
    @eval if hasmethod(∇Φ,($f,Real)) # check if this method is defined for the current type
        @eval ∇Φ(uu::u.AccelerationUnits,$f,args...;kws...) = ∇Φ($f,args...;kws...) * defaultunits.∇Φunit |> uu
        @eval ∇Φ(uu::u.AccelerationUnits,$f,r::u.Length,args...;kws...) = ∇Φ($f,homogenize_units(r),args...;kws...) * defaultunits.∇Φunit |> uu
        @eval ∇Φ($f,r::u.Length,args...;kws...) = ∇Φ($f,homogenize_units(r),args...;kws...) * defaultunits.∇Φunit
    end
    @eval if hasmethod(∇∇Φ,($f,Real)) # check if this method is defined for the current type
        @eval ∇∇Φ(uu::∇∇ΦdimensionUnits,$f,args...;kws...) = ∇∇Φ($f,args...;kws...) * defaultunits.∇∇Φunit |> uu
        @eval ∇∇Φ(uu::∇∇ΦdimensionUnits,$f,r::u.Length,args...;kws...) = ∇∇Φ($f,homogenize_units(r),args...;kws...) * defaultunits.∇∇Φunit |> uu
        @eval ∇∇Φ($f,r::u.Length,args...;kws...) = ∇∇Φ($f,homogenize_units(r),args...;kws...) * defaultunits.∇∇Φunit
    end
end
