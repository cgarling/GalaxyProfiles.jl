import .UnitfulAstro as ua
import .Unitful as u

# Define a dimensionality for surface densities
u.@derived_dimension SurfaceDensity u.𝐌/u.𝐋^2
u.@derived_dimension ∇ρdimension u.𝐌/u.𝐋^4
u.@derived_dimension Φdimension u.𝐋^2/u.𝐓^2
u.@derived_dimension ∇mdimension u.𝐌/u.𝐋
# u.@derived_dimension ∇Φdimension u.𝐋/u.𝐓^2  # this is just u.AccelerationUnits

# module to hold default units
module defaultunits
import ..UnitfulAstro as ua
import ..Unitful as u

const mass = ua.Msun
const ∇mass = ua.Msun / ua.kpc
const density = ua.Msun / ua.kpc^3
const ∇density = ua.Msun / ua.kpc^4
const surfacedensity = ua.Msun / ua.kpc^2
const length = ua.kpc
const velocity = u.km / u.s
const Φunit = u.km^2 / u.s^2
const ∇Φunit = u.km^2 / u.s^2 / ua.kpc
const ∇∇Φunit = u.km^2 / u.s^2 / ua.kpc^2
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
ExponentialDisk(Σ0::SurfaceDensity,rs::u.Length) = ExponentialDisk(homogenize_units(Σ0),homogenize_units(rs))
function ExponentialDisk(rs::u.Length;M=nothing,Σ0=nothing)
    if isnothing(M)
        @assert !isnothing(Σ0)
        @assert Σ0 isa SurfaceDensity
        ExponentialDisk(Σ0,rs)
    else
        @assert M isa u.Mass
        ExponentialDisk(M/(2π*rs^2),rs)
    end
end

GeneralIsothermal(ρ0::u.Density,rs::u.Length,α::Real) = GeneralIsothermal(homogenize_units(ρ0),homogenize_units(rs),α)
GeneralIsothermal(rs::u.Length,α::Real,M::u.Mass,Rmax::u.Length) = GeneralIsothermal(homogenize_units(rs),α,homogenize_units(M),homogenize_units(Rmax))

SIS(ρ0::u.Density,rs::u.Length) = SIS(homogenize_units(ρ0),homogenize_units(rs))
SIS(rs::u.Length,M::u.Mass,Rmax::u.Length) = SIS(homogenize_units(rs),homogenize_units(M),homogenize_units(Rmax))

NFW(ρ0::u.Density,rs::u.Length) = NFW(homogenize_units(ρ0),homogenize_units(rs))

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
for f in (:ExponentialDisk,:GeneralIsothermal) # common quantities for both 3D densities and 2D SurfaceDensities
    @eval Σ(uu::SurfaceDensityUnits,$f,args...;kws...) = Σ($f,args...;kws...) * defaultunits.surfacedensity |> uu
    @eval Σ(uu::SurfaceDensityUnits,$f,r::u.Length,args...;kws...) = Σ($f,homogenize_units(r),args...;kws...) * defaultunits.surfacedensity |> uu
    @eval Σ($f,r::u.Length,args...;kws...) = Σ($f,homogenize_units(r),args...;kws...) * defaultunits.surfacedensity
    @eval ∇Σ(uu::u.DensityUnits,$f,args...;kws...) = ∇Σ($f,args...;kws...) * defaultunits.density |> uu
    @eval ∇Σ(uu::u.DensityUnits,$f,r::u.Length,args...;kws...) = ∇Σ($f,homogenize_units(r),args...;kws...) * defaultunits.density |> uu
    @eval ∇Σ($f,r::u.LengthUnits,args...;kws...) = ∇Σ($f,homogenize_units(r),args...;kws...) * defaultunits.density
    # Σmean
    @eval invΣ(uu::u.LengthUnits,$f,args...;kws...) = invΣ($f,args...;kws...) * defaultunits.length |> uu
    @eval invΣ(uu::u.LengthUnits,$f,x::SurfaceDensity,args...;kws...) = invΣ($f,homogenize_units(x),args...;kws...) * defaultunits.length |> uu
    @eval invΣ($f,x::SurfaceDensity,args...;kws...) = invΣ($f,homogenize_units(x),args...;kws...) * defaultunits.length
end
for f in (:GeneralIsothermal,) # quantities requiring 3D densities
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
    # ρmean
    # invρmean
    # Σ is above
    # ∇Σ is above
    # Σmean is above
    # invΣ is above
    @eval M(uu::u.MassUnits,$f,args...;kws...) = M($f,args...;kws...) * defaultunits.mass |> uu
    @eval M(uu::u.MassUnits,$f,r::u.Length,args...;kws...) = M($f,homogenize_units(r),args...;kws...) * defaultunits.mass |> uu
    @eval M($f,r::u.Length,args...;kws...) = M($f,homogenize_units(r),args...;kws...) * defaultunits.mass
    @eval ∇M(uu::∇mdimensionUnits,$f,args...;kws...) = ∇M($f,args...;kws...) * defaultunits.∇mass |> uu
    @eval ∇M(uu::∇mdimensionUnits,$f,r::u.Length,args...;kws...) = ∇M($f,homogenize_units(r),args...;kws...) * defaultunits.∇mass |> uu
    @eval ∇M($f,r::u.Length,args...;kws...) = ∇M($f,homogenize_units(r),args...;kws...) * defaultunits.∇mass
    @eval invM(uu::u.LengthUnits,$f,args...;kws...) = invM($f,args...;kws...) * defaultunits.length |> uu
    @eval invM(uu::u.LengthUnits,$f,x::u.Mass,args...;kws...) = invM($f,homogenize_units(x),args...;kws...) * defaultunits.length |> uu
    @eval invM($f,x::u.Mass,args...;kws...) = invM($f,homogenize_units(x),args...;kws...) * defaultunits.length
    # Mtot
    @eval Φ(uu::ΦdimensionUnits,$f,args...;kws...) = Φ($f,args...;kws...) * defaultunits.Φunit |> uu
    @eval Φ(uu::ΦdimensionUnits,$f,r::u.Length,args...;kws...) = Φ($f,homogenize_units(r),args...;kws...) * defaultunits.Φunit |> uu
    @eval Φ($f,r::u.Length,args...;kws...) = Φ($f,homogenize_units(r),args...;kws...) * defaultunits.Φunit
    @eval ∇Φ(uu::u.AccelerationUnits,$f,args...;kws...) = ∇Φ($f,args...;kws...) * defaultunits.∇Φunit |> uu
    @eval ∇Φ(uu::u.AccelerationUnits,$f,r::u.Length,args...;kws...) = ∇Φ($f,homogenize_units(r),args...;kws...) * defaultunits.∇Φunit |> uu
    @eval ∇Φ($f,r::u.Length,args...;kws...) = ∇Φ($f,homogenize_units(r),args...;kws...) * defaultunits.∇Φunit
end
# this ALMOST works ... need to do functions with different dimensionalities separately, or define some
# functions that will return the correct units depending on the function or the output
# for t in (:ExponentialDisk,:GeneralIsothermal)
#     for f in (:Σ,:∇Σ)
#         @eval $f(uu::u.Unitlike,$t,args...;kws...) = $f($t,args...;kws...) * u.ustrip(uu,1*ua.Msun/ua.kpc^2)
#         @eval $f(uu::u.Unitlike,$t,r::u.Length,args...;kws...) = $f($t,u.ustrip(ua.kpc,r),args...;kws...) * u.ustrip(uu,1*ua.Msun/ua.kpc^2)
#         @eval $f($t,r::u.Length,args...;kws...) = $f($t,u.ustrip(ua.kpc,r),args...;kws...)
#     end
# end
