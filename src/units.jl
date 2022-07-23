import .UnitfulAstro as ua
import .Unitful as u

# Define a dimensionality for surface densities
u.@derived_dimension SurfaceDensity u.𝐌/u.𝐋^2

# Dictionary to hold default units
# this is slow
# const default_units=Dict("length"=>ua.kpc)

# General functions for converting units
# Should check out how to have user-controllable default units, rather than the hard-coded values here
convert_density(ρ0::u.Density) = u.ustrip(ua.Msun/ua.kpc^3, ρ0)
convert_mass(M::u.Mass) = u.ustrip(ua.Msun, M)
convert_surface_density(Σ0::SurfaceDensity) = u.ustrip(ua.Msun/ua.kpc^2, Σ0)
convert_length(r::u.Length) = u.ustrip(ua.kpc, r)
# convert_length(r::u.Length) = u.ustrip(r|>default_units["length"])

# Create Unitful constructors for our various composite types
ExponentialDisk(Σ0::SurfaceDensity,rs::u.Length) = ExponentialDisk(convert_surface_density(Σ0),convert_length(rs))
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

GeneralIsothermal(ρ0::u.Density,rs::u.Length,α::Real) = GeneralIsothermal(convert_density(ρ0),convert_length(rs),α)
GeneralIsothermal(rs::u.Length,α::Real,M::u.Mass,Rmax::u.Length) = GeneralIsothermal(convert_length(rs),α,convert_mass(M),convert_length(Rmax))

SIS(ρ0::u.Density,rs::u.Length) = SIS(convert_density(ρ0),convert_length(rs))
SIS(rs::u.Length,M::u.Mass,Rmax::u.Length) = SIS(convert_length(rs),convert_mass(M),convert_length(Rmax))

NFW(ρ0::u.Density,rs::u.Length) = NFW(convert_density(ρ0),convert_length(rs))

#########################################################################################
# Easily select a different unit
# This is pretty annoying but it's nice to have I guess ...
# This also has hard-coded default units, which I should probably specify somewhere explicitly
# Switch this to have loops over Types and over Functions

for f in (:ExponentialDisk,:GeneralIsothermal) # common quantities for both 3D densities and 2D SurfaceDensities
    @eval Σ(uu::u.Unitlike,$f,args...;kws...) = Σ($f,args...;kws...) * u.ustrip(1*ua.Msun/ua.kpc^2 |> uu)
    @eval Σ(uu::u.Unitlike,$f,r::u.Length,args...;kws...) = Σ($f,u.ustrip(r|>ua.kpc),args...;kws...) * u.ustrip(1*ua.Msun/ua.kpc^2 |> uu)
    @eval Σ($f,r::u.Length,args...;kws...) = Σ($f,u.ustrip(r|>ua.kpc),args...;kws...)
    @eval ∇Σ(uu::u.Unitlike,$f,args...;kws...) = ∇Σ($f,args...;kws...) * u.ustrip(1*ua.Msun/ua.kpc^3 |> uu)
    @eval ∇Σ(uu::u.Unitlike,$f,r::u.Length,args...;kws...) = ∇Σ($f,u.ustrip(r|>ua.kpc),args...;kws...) * u.ustrip(1*ua.Msun/ua.kpc^3 |> uu)
    @eval ∇Σ($f,r::u.Length,args...;kws...) = ∇Σ($f,u.ustrip(r|>ua.kpc),args...;kws...)
end
for f in (:GeneralIsothermal,) # quantities requiring 3D densities
    @eval Φ(uu::u.Unitlike,$f,args...;kws...) = Φ($f,args...;kws...) * u.ustrip(1*u.km^2/u.s^2 |> uu)
    @eval Φ(uu::u.Unitlike,$f,r::u.Length,args...;kws...) = Φ($f,u.ustrip(r|>ua.kpc),args...;kws...) * u.ustrip(1*u.km^2/u.s^2 |> uu)
    @eval Φ($f,r::u.Length,args...;kws...) = Φ($f,u.ustrip(r|>ua.kpc),args...;kws...)
    @eval ∇Φ(uu::u.Unitlike,$f,args...;kws...) = ∇Φ($f,args...;kws...) * u.ustrip(1*u.km^2/u.s^2/ua.kpc |> uu)
    @eval ∇Φ(uu::u.Unitlike,$f,r::u.Length,args...;kws...) = ∇Φ($f,u.ustrip(r|>ua.kpc),args...;kws...) * u.ustrip(1*u.km^2/u.s^2/ua.kpc |> uu)
    @eval ∇Φ($f,r::u.Length,args...;kws...) = ∇Φ($f,u.ustrip(r|>ua.kpc),args...;kws...)
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
