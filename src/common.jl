module constants
# const Gkpc = 4.5171030498949647e-39     # G in kpc^3/solMass/s^2
# const Gpc = 4.5171030498949647e-30      # G in pc^3/solMass/s^2
# const Gvelkpc = 4.30091727003628e-6     # G in kpc*km^2/solMass/s^2, for velocity calculations
# const Gvelkpc2 = 1.3938323614347172e-22 # G in kpc^2*km/solMass/s^2, for velocity calculations
# const Gvelpc = 4.30091727003628e-3      # G in kpc*km^2/solMass/s^2
# Use BigFloat for more accuracy
const Gkpc = big"4.517103049894965029489870560544420380775973240960280169677794680278054678900598e-39" # G in kpc^3/solMass/s^2
const Gpc = big"4.517103049894965632856118045698970338735168159241032233763499248446748879359808e-30"  # G in pc^3/solMass/s^2
const Gvelkpc = big"4.3009172700362805113763897679746150970458984375e-6"   # G in kpc*km^2/solMass/s^2, for velocity and Φ calculations; (G |> ua.kpc*u.km^2/ua.Msun/u.s^2)
const Gvelkpc2 = big"1.393832361434717359905559114367183944790625806098294248158708796836435794830327e-22" # G in kpc^2*km/solMass/s^2, for ∇Φ calculations; (G |> ua.kpc^2*u.km/ua.Msun/u.s^2)
const Gvelpc = big"4.3009172700362805113763897679746150970458984375e-6"
const nfwvmaxpar = big"2.162581587064609834856553669603264573507154023820515448359231774428038693093796" # constant for NFW's Vmax
end # module

module utilities
""" `get_inf(x::T) where T` returns the infinity of type `T` if `T` is `Float64`, `Float32`, or `Float16`. Otherwise, returns `Inf`."""
get_inf(::Float64) = Inf64
get_inf(::Float32) = Inf32
get_inf(::Float16) = Inf16
get_inf(::Any) = Inf
# Useful but not currently used
# function hasspecialization(f::F, d::T, r::S) where {F, T, S}
#     m = which(f, Tuple{T, S}) # 
#     return m.sig == Tuple{F, T.name.wrapper, Real}
# end
end # module

"""
    params(d::AbstractMassProfile)

Get the fields of the struct as expected by its defined methods.
"""
params(d::AbstractMassProfile)
"""
    scale_radius(d::AbstractMassProfile)
    scale_radius(uu::Unitful.LengthUnits, d::AbstractMassProfile)

Returns the characteristic scale radius of the profile; used for some default methods. An example is `rs` for the [`ExponentialDisk`](@ref) model. 
"""
scale_radius(d::AbstractMassProfile)
"""
    ρ(d::AbstractDensity, r::Real)
    ρ(uu::Unitful.DensityUnits, d::AbstractDensity, r::Real)
    ρ(d::AbstractDensity, r::Unitful.Length)
    ρ(uu::Unitful.DensityUnits, d::AbstractDensity, r::Unitful.Length)

Returns the density of profile `d` [M⊙ kpc^-3] at radius `r` [kpc].
"""
ρ(d::AbstractDensity, r::Real) # ``\\left[ \\text{M}_\\odot \\, \\text{kpc}^{-3} \\right]``
"""
    invρ(d::AbstractDensity, x::Real) # Only valid if specialized method exists
    invρ(d::AbstractDensity, x::Real,
         interval::NTuple{2,Real}=(scale_radius(d)/100,
                                   100*scale_radius(d)); kws...) # Root-finding generic method
    invρ(uu::Unitful.LengthUnits, d::AbstractDensity, x::Real)
    invρ(d::AbstractDensity, x::Unitful.Density)
    invρ(uu::Unitful.LengthUnits, d::AbstractDensity, x::Unitful.Density)

Solve for the radius `r` [kpc] at which the density [M⊙ kpc^-3] is `x` for density profile `d`. Requires `x>0`.

When this method is not specialized for `d`, it will use an interval bracketing method from [`Roots.jl`](https://github.com/JuliaMath/Roots.jl), requiring that `ρ(d,r)` be defined. For this method, `kws...` are passed to `Roots.find_zero`.
"""
function invρ(d::AbstractDensity{T}, x::S,
              interval::NTuple{2,V} =
                  (scale_radius(d)/100, 100*scale_radius(d));
              kws...) where {T <: Real, S <: Real, V <: Real}
    @assert x > 0 "x must be > 0"
    U = promote_type(T, S, V)
    return find_zero(y->ρ(d,y)-x, (U(interval[1]), U(interval[2])); kws...)
end
"""
    ∇ρ(d::AbstractDensity, r::Real)
    ∇ρ(uu::GalaxyProfiles.∇ρdimensionUnits, d::AbstractDensity, r::Real)
    ∇ρ(d::AbstractDensity, r::Unitful.Length)
    ∇ρ(uu::GalaxyProfiles.∇ρdimensionUnits, d::AbstractDensity, r::Unitful.Length)

Returns the radial derivative of the density [M⊙ kpc^-4] of the mass profile `d` evaluated at radius `r` [kpc].
"""
∇ρ(d::AbstractDensity, r::Real)
"""
    ρmean(d::AbstractDensity, r::Real)
    ρmean(uu::Unitful.DensityUnits, d::AbstractDensity, r::Real)
    ρmean(d::AbstractDensity, r::Unitful.Length)
    ρmean(uu::Unitful.DensityUnits, d::AbstractDensity, r::Unitful.Length)

Returns the average density [M⊙ kpc^-3] of density profile `d` inside `r` [kpc]. The generic method uses enclosed mass divided by spherical volume; `M(d,r) * 3 / (4π*r^3)`. 
"""
ρmean(d::AbstractDensity, r::Real) = M(d,r) / r^3 * 3 / 4 / π
"""
    invρmean(d::AbstractDensity, x::Real)
    invρmean(d::AbstractDensity, x::Real,
             interval::NTuple{2,Real}=(scale_radius(d)/100,
                                       100*scale_radius(d)); kws...)
    invρmean(uu::Unitful.LengthUnits, d::AbstractDensity, x::Real)
    invρmean(d::AbstractDensity, x::Unitful.Density)
    invρmean(uu::Unitful.LengthUnits, d::AbstractDensity, x::Unitful.Density)

Solve for the radius `r` [kpc] inside which the average density is `x` [M⊙ kpc^-3]. Requires `x>0`.

When this method is not specialized for `d`, it will use an interval bracketing method from [`Roots.jl`](https://github.com/JuliaMath/Roots.jl), requiring that `ρmean(d,r)` or `M(d,r)` be defined. For this method, `kws...` are passed to `Roots.find_zero`.
"""
function invρmean(d::AbstractDensity{T}, x::S,
                  interval::NTuple{2,V} =
                      (scale_radius(d)/100, 100*scale_radius(d));
                  kws...) where {T <: Real, S <: Real, V <: Real}
    @assert x > 0 "x must be > 0"
    U = promote_type(T, S, V)
    return find_zero(y->ρmean(d,y)-x, (U(interval[1]), U(interval[2])); kws...)
end
"""
    Σ(d::AbstractMassProfile, r::Real)
    Σ(uu::GalaxyProfiles.SurfaceDensityUnits, d::AbstractMassProfile, r::Real)
    Σ(d::AbstractMassProfile, r::Unitful.Length)
    Σ(uu::GalaxyProfiles.SurfaceDensityUnits, d::AbstractMassProfile, r::Unitful.Length)

Returns the surface density [M⊙ kpc^-2] of mass profile `d` at radius `r`. For 3D density profiles (i.e., subtypes of [`AbstractDensity`](@ref GalaxyProfiles.AbstractDensity)), this will be the *projected* surface density, which, for spherical systems, is defined by the Abel integral 

```math
\\Sigma(r) = 2 \\int_R^\\infty \\rho(r) \\frac{r}{ \\sqrt{r^2-R^2} } \\, dr
```

which has an inverse of

```math
\\rho(r) = -\\frac{1}{\\pi} \\int_r^\\infty \\frac{d\\Sigma(R)}{dR} \\frac{dR}{\\sqrt{R^2-r^2}}
```
"""
Σ(d::AbstractMassProfile, r::Real) = 2 * quadgk(x->x*ρ(d,x)/sqrt(x^2-r^2), r, utilities.get_inf(r))[1]
"""
    ∇Σ(d::AbstractMassProfile, r::Real)
    ∇Σ(uu::Unitful.DensityUnits, d::AbstractMassProfile, r::Real)
    ∇Σ(d::AbstractMassProfile, r::Unitful.Length)
    ∇Σ(uu::Unitful.DensityUnits, d::AbstractMassProfile, r::Unitful.Length)

Returns the radial derivative of the surface density [M⊙ kpc^-3] of the mass profile `d` evaluated at radius `r` [kpc].
"""
∇Σ(d::AbstractMassProfile, r::Real)
"""
    Σmean(d::AbstractMassProfile, r::Real)
    Σmean(uu::GalaxyProfiles.SurfaceDensityUnits, d::AbstractMassProfile, r::Real)
    Σmean(d::AbstractMassProfile, r::Unitful.Length)
    Σmean(uu::GalaxyProfiles.SurfaceDensityUnits, d::AbstractMassProfile, r::Unitful.Length)

Returns the mean projected surface density [M⊙ kpc^-2] of the mass profile `d` inside the radius `r` [kpc]. The generic method uses projected enclosed mass divided by area: `Mproj(d::AbstractMassProfile,r::Real) / (π * r^2)`.
"""
Σmean(d::AbstractMassProfile, r::Real) = Mproj(d,r) / r^2 / π
"""
    invΣ(d::AbstractMassProfile, x::Real)
    invΣ(d::AbstractMassProfile, x::Real,
         interval::NTuple{2,Real}=(scale_radius(d)/100,
                                   100*scale_radius(d)); kws...)
    invΣ(uu::Unitful.LengthUnits, d::AbstractMassProfile, x::Real)
    invΣ(d::AbstractMassProfile, r::Unitful.Length)
    invΣ(uu::Unitful.LengthUnits, d::AbstractMassProfile, r::Unitful.Length)

Solve for the radius `r` [kpc] at which the surface density [M⊙ kpc^-2] is `x` for profile `d`. Requires `x>0`.

When this method is not specialized for `d`, it will use an interval bracketing method from [`Roots.jl`](https://github.com/JuliaMath/Roots.jl), requiring that `Σ(d,r)` be defined. For this method, `kws...` are passed to `Roots.find_zero`.
"""
function invΣ(d::AbstractMassProfile{T}, x::S,
              interval::NTuple{2,V} =
                  (scale_radius(d)/100, 100*scale_radius(d));
              kws...) where {T <: Real, S <: Real, V <: Real}
    @assert x > 0 "x must be > 0"
    U = promote_type(T, S, V)
    return find_zero(y->Σ(d,y)-x, (U(interval[1]), U(interval[2])); kws...)
end
"""
    M(d::AbstractDensity, r::Real)
    M(uu::Unitful.MassUnits, d::AbstractDensity, r::Real)
    M(d::AbstractDensity, r::Unitful.Length)
    M(uu::Unitful.MassUnits, d::AbstractDensity, r::Unitful.Length)

Returns the total mass [M⊙] of the profile `d` enclosed within a spherical shell of radius `r`. For spherical systems this is given by the integral

```math
M(\\lt R) = 4\\pi \\int_0^R r^2 ρ(r) dr
```

When this method is not specialized for `d`, it will compute the numerical integral using `QuadGK.quadgk`.
"""
M(d::AbstractDensity, r::Real) = quadgk(x->x^2 * ρ(d,x), zero(r), r)[1] * 4 * π
"""
    ∇M(d::AbstractDensity, r::Real)
    ∇M(uu::GalaxyProfiles.∇mdimensionUnits, d::AbstractDensity, r::Real)
    ∇M(d::AbstractDensity, r::Unitful.Length)
    ∇M(uu::GalaxyProfiles.∇mdimensionUnits, d::AbstractDensity, r::Unitful.Length)

The radial derivative of the enclosed mass [M⊙ kpc^-1] of profile `d` evaluated at radius `r`.
"""
∇M(d::AbstractDensity, r::Real)
"""
    invM(d::AbstractDensity, x::Real)
    invM(d::AbstractDensity, x::Real,
         interval::NTuple{2,Real}=(scale_radius(d)/100,
                                   100*scale_radius(d)); kws...)
    invM(uu::Unitful.LengthUnits, d::AbstractDensity, x::Real)
    invM(d::AbstractDensity, x::Unitful.Mass)
    invM(uu::Unitful.LengthUnits, d::AbstractDensity, x::Unitful.Mass)

Solve for the radius `r` [kpc] at which the enclosed mass [M⊙] is `x` for profile `d`. Requires `x>0`.

When this method is not specialized for `d`, it will use an interval bracketing method from [`Roots.jl`](https://github.com/JuliaMath/Roots.jl), requiring that `M(d,r)` be defined. For this method, `kws...` are passed to `Roots.find_zero`.
"""
function invM(d::AbstractDensity{T}, x::S,
              interval::NTuple{2,V} =
                  (scale_radius(d)/100, 100*scale_radius(d));
              kws...) where {T <: Real, S <: Real, V <: Real}
    @assert x > 0 "x must be > 0"
    U = promote_type(T, S, V)
    return find_zero(y->M(d,y)-x, (U(interval[1]), U(interval[2])); kws...)
end
"""
    Mtot(d::AbstractMassProfile)
    Mtot(uu::Unitful.MassUnits, d::AbstractMassProfile)

Returns the total mass [kpc] of the profile `d`. For profiles which do not have a well-defined total mass `M(d,Inf)`, this quantity should be defined as the limit ``\\lim_{r \\to \\infty} M(r)``.
"""
Mtot(d::AbstractMassProfile)
"""
    Mproj(d::AbstractMassProfile, r::Real)
    Mproj(uu::Unitful.MassUnits, d::AbstractMassProfile, r::Real)
    Mproj(d::AbstractMassProfile, r::Unitful.Length)
    Mproj(uu::Unitful.MassUnits, d::AbstractMassProfile, r::Unitful.Length)

Returns the total line-of-sight projected mass [M⊙] enclosed within a radius `r` [kpc] for the profile `d`. For spherical systems this is given by the integral

```math
M_{\\text{proj}}(\\lt R) = 2\\pi \\int_0^R r \\, \\Sigma(r) dr = 2\\pi \\int_0^R r \\left( 2 \\int_r^\\infty \\rho(r^\\prime) \\frac{r^\\prime}{ \\sqrt{r^{\\prime \\, 2}-r^2} } \\, dr^\\prime \\right) dr
```
"""
Mproj(d::AbstractMassProfile, r::Real) = quadgk(x->x * Σ(d,x), zero(r), r)[1] * 2 * π
"""
    ∇Mproj(d::AbstractMassProfile, r::Real)
    ∇Mproj(uu::GalaxyProfiles.∇mdimensionUnits, d::AbstractMassProfile, r::Real)
    ∇Mproj(d::AbstractMassProfile, r::Unitful.Length)
    ∇Mproj(uu::GalaxyProfiles.∇mdimensionUnits, d::AbstractMassProfile, r::Unitful.Length)

The radial derivative of the line-of-sight projected mass [M⊙ kpc^-1] of profile `d` enclosed within a radius `r` [kpc].
"""
∇Mproj(d::AbstractMassProfile, r::Real)
"""
    invMproj(d::AbstractMassProfile, x::Real)
    invMproj(d::AbstractMassProfile, x::Real,
             interval::NTuple{2,Real}=(scale_radius(d)/100,
                                       100*scale_radius(d)); kws...)
    invMproj(d::AbstractMassProfile, x::Unitful.Mass)
    invMproj(uu::Unitful.LengthUnits, d::AbstractMassProfile, x::Unitful.Mass)

Solve for the radius `r` [kpc] at which the line-of-sight projected enclosed mass is `x` [M⊙] for profile `d`. Requires `x>0`.

When this method is not specialized for `d`, it will use an interval bracketing method from [`Roots.jl`](https://github.com/JuliaMath/Roots.jl), requiring that `M(d,r)` be defined. For this method, `kws...` are passed to `Roots.find_zero`.
"""
function invMproj(d::AbstractMassProfile{T}, x::S,
                  interval::NTuple{2,V} =
                      (scale_radius(d)/100, 100*scale_radius(d));
                  kws...) where {T <: Real, S <: Real, V <: Real}
    @assert x > 0 "x must be > 0"
    U = promote_type(T, S, V)
    return find_zero(y->Mproj(d,y)-x, (U(interval[1]), U(interval[2])); kws...)
end
"""
    dynamical_time(d::Union{AbstractDensity, AbstractMassProfile}, r::Real)
    dynamical_time(d::Union{AbstractDensity, AbstractMassProfile}, r::Unitful.Length)
    dynamical_time(uu::Unitful.TimeUnits, d::Union{AbstractDensity, AbstractMassProfile}, r::Unitful.Length)

Returns the dynamical time [yr] at radius `r` [kpc] in the provided density profile `d`. The dynamical time is a characteristic timescale associated with orbits in potentials. This implementation specifically calculates Equation 2.39 in Binney & Tremaine Galactic Dynamics 2E, with the replacement of the standard density [`ρ`](@ref) for the average density interior to `r`, [`ρmean`](@ref), to better account for inhomogenous systems. 

```math
t_\\text{dyn} \\left( r \\right) = \\sqrt{ \\frac{3\\pi}{16 \\, G \\, \\overline{\\rho} \\left( r \\right)} }
```

The above equation is used when the argument `d::AbstractDensity`. As `ρmean` can be expressed as a function of the total mass interior to the orbit, this can also be written as

```math
\\begin{aligned}
\\overline{\\rho} \\left( r \\right) &= \\frac{3 \\, M \\left( r \\right)}{4 \\, \\pi \\, r^3} \\newline
t_\\text{dyn} \\left( r \\right) &= \\pi \\, \\sqrt{ \\frac{r^3}{4 \\, G \\, M \\left( r \\right)} }
\\end{aligned}
```

where ``M \\left( r \\right)`` is the mass enclosed inside radius ``r``, provided by the method [`M`](@ref). This equation is used when the argument `d::AbstractMassProfile`, with the substitution of the projected mass [`Mproj`](@ref) for the 3-D enclosed mass [`M`](@ref), which is undefined for profiles that do not have `ρ` densities. 

# Alternative Definitions
Note that some authors prefer to omit the factor of ``\\sqrt{ \\frac{1}{16} } = 1/4`` in the denominator of the first equation above; this is, for example, the convention used by galpy. Dynamical times thus defined will be four times larger than those calculated by this method. 
"""
function dynamical_time(d::AbstractDensity{T}, r::S) where {T <: Real, S <: Real}
    U = promote_type(T, S)
    # constants.Gkpc has time units of 1/s^2; 31557600 is seconds in year to return value in yr
    # return sqrt(inv(GalaxyProfiles.ρmean(d, r)) * 3 * π / 16 / T(constants.Gkpc)) / 31557600
    return sqrt(inv(GalaxyProfiles.ρmean(d, r)) * 3 * π) / (4 * 31557600 * sqrt(U(constants.Gkpc)))
end
dynamical_time(d::AbstractMassProfile{T}, r::S) where {T <: Real, S <: Real} =
    π * sqrt(r^3 / Mproj(d, r)) / (2 * 31557600 * sqrt(convert(promote_type(T,S), constants.Gkpc)))
    # π * sqrt( r^3 / 4 / T(constants.Gkpc) / M(d, r) ) / 31557600
"""
    cdf2D(d::AbstractMassProfile, r::Real)
    cdf2D(d::AbstractMassProfile, r::Unitful.Quantity)

Returns the value of the cumulative distribution function of the profile `d` at `r` [kpc] in two dimensions (i.e., along a line of sight). This is defined as `Σ(d,r) / Mtot(d)`.
"""
cdf2D(d::AbstractMassProfile, r::Real) = Σ(d,r) / Mtot(d)
"""
    cdf3D(d::AbstractDensity, r::Real)
    cdf3D(d::AbstractDensity, r::Unitful.Quantity)

Returns the value of the cumulative distribution function of the profile `d` at `r` [kpc] in three dimensions. This is defined as `M(d,r) / Mtot(d)`.
"""
cdf3D(d::AbstractDensity, r::Real) = M(d,r) / Mtot(d)
"""
    ccdf2D(d::AbstractMassProfile, r::Real)
    ccdf2D(d::AbstractMassProfile, r::Unitful.Quantity)

Returns the value of the complementary cumulative distribution function of the profile `d` at `r` [kpc] in two dimensions (i.e., along a line of sight). This is defined as `1 - cdf2D(d,r) = 1 - Σ(d,r) / Mtot(d)`. 
"""
ccdf2D(d::AbstractMassProfile, r::Real) = 1 - cdf2D(d,r)
"""
    ccdf3D(d::AbstractMassProfile, r::Real)
    ccdf3D(d::AbstractMassProfile, r::Unitful.Quantity)

Returns the value of the complementary cumulative distribution function of the profile `d` at `r` [kpc] in three dimensions. This is defined as `1 - cdf3D(d,r) = 1 - M(d,r) / Mtot(d)`. 
"""
ccdf3D(d::AbstractMassProfile, r::Real) = 1 - cdf3D(d,r)
"""
    quantile2D(d::AbstractMassProfile, r::Real)
    quantile2D(d::AbstractMassProfile, r::Unitful.Quantity)

Returns the value of the inverse cumulative distribution function of the profile `d` at `r` [kpc] in two dimensions (i.e., along a line of sight).
"""
quantile2D(d::AbstractMassProfile, r::Real)
"""
    quantile3D(d::AbstractDensity, r::Real)
    quantile3D(d::AbstractDensity, r::Unitful.Quantity)

Returns the value of the inverse cumulative distribution function of the profile `d` at `r` [kpc] in three dimensions.
"""
quantile3D(d::AbstractDensity, r::Real)
"""
    cquantile2D(d::AbstractMassProfile, r::Real)
    cquantile2D(d::AbstractMassProfile, r::Unitful.Quantity)

Returns the value of the complementary quantile of the profile `d` at `r` [kpc] in two dimensions (i.e., along a line of sight). This is defined as `quantile2D(d, 1-r)`.
"""
cquantile2D(d::AbstractMassProfile, r::Real)
"""
    cquantile3D(d::AbstractDensity, r::Real)
    cquantile3D(d::AbstractDensity, r::Unitful.Quantity)

Returns the value of the complementary quantile of the profile `d` at `r` [kpc] in three dimensions. This is defined as `quantile3D(d, 1-r)`.
"""
cquantile3D(d::AbstractDensity, r::Real)
"""
    Vcirc(d::AbstractDensity, r::Real)
    Vcirc(uu::Unitful.VelocityUnits, d::AbstractDensity, r::Real)
    Vcirc(d::AbstractDensity, r::Unitful.Length)
    Vcirc(uu::Unitful.VelocityUnits, d::AbstractDensity, r::Unitful.Length)

Returns the circular velocity [km/s] at `r` [kpc], defined as the speed of a particle of insignificant mass in a circular orbit at radius `r`. This is calculated as

```math
v_c^2(r) = \\frac{G M(r)}{r} = r \\frac{d\\Phi}{dr} = r\\nabla\\Phi(r)
```

By default uses `G` in units such that if `rs` and `r` are in kpc, the velocity ends up in `km/s`. For example, for [`GeneralIsothermal`](@ref) we have `[G] = [kpc * km^2 / Msun / s^2]` so that the velocity ends up in `km/s`. This method has a generic implementation of `sqrt( GalaxyProfiles.constants.Gvelkpc * M(d::AbstractDensity,r) / r)`.
"""
Vcirc(d::AbstractDensity{T}, r::S) where {T <: Real, S <: Real} =
    sqrt(convert(promote_type(T,S), constants.Gvelkpc) * M(d, r) / r)
"""
    Vesc(d::AbstractDensity, r::Real)
    Vesc(uu::Unitful.VelocityUnits, d::AbstractDensity, r::Real)
    Vesc(d::AbstractDensity, r::Unitful.Length)
    Vesc(uu::Unitful.VelocityUnits, d::AbstractDensity, r::Unitful.Length)

Returns the escape velocity [km/s] at `r` [kpc], calculated as
```math
v^2_{\\text{esc}}(r) = 2 |\\Phi(r)|
```
if ``\\Phi \\to 0`` for ``r \\to \\infty``; see the note for [`Φ`](@ref).

By default uses `G` in units such that if `rs` and `r` are in kpc, the velocity ends up in `km/s`. For example, for [`GeneralIsothermal`](@ref) we have `[G] = [kpc * km^2 / Msun / s^2]` so that the velocity ends up in `km/s`. 
"""
Vesc(d::AbstractDensity, r::Real) = sqrt(2 * abs(Φ(d, r)))
"""
    Vmax(d::AbstractDensity)
    Vmax(uu::Unitful.VelocityUnits, d::AbstractDensity)

Returns the maximum circular velocity [km/s] of `d` and the corresponding radius [kpc]. Can be found by solving

```math
    \\frac{d v_c(r)}{dr} = 0
```

for `r`, where ``v_c`` is the circular velocity, then evaluating the circular velocity at `r`.
"""
Vmax(d::AbstractDensity)
"""
    σr(d::AbstractDensity, r::Real, β::Real)

Returns the radial velocity dispersion [km/s] of `d` at radius `r` [kpc] for constant velocity anisotropy `β` given by Equation 4.216 in Binney & Tremaine Galactic Dynamics 2E,

```math
\\sigma_r^2 \\left( R \\right) = \\frac{1}{R^{2\\beta} \\, \\rho(R)} \\int_R^\\infty r^{2\\beta} \\, \\rho(r) \\, \\frac{d\\Phi}{dr} \\, dr
```

and as ``\\frac{d\\Phi}{dr} = -G M(r) / r^2`` we can alternatively write

```math
\\sigma_r^2 \\left( R \\right) = \\frac{G}{R^{2\\beta} \\, \\rho(R)} \\int_R^\\infty r^{2\\left( \\beta-1 \\right)} \\, \\rho(r) \\, M\\left( r \\right) \\, dr
```

which is the generic method as we expect ``M(r)`` to be more commonly available than ``\\frac{d\\Phi}{dr}``.
"""
σr(d::AbstractDensity, r::T, β::S) where {T <: Real, S <: Real} = σr(d, promote(r, β)...)
function σr(d::AbstractDensity{T}, r::S, β::S) where {T <: Real, S <: Real}
    U = promote_type(T, S)
    return sqrt(U(constants.Gvelkpc) / r^(2β) / ρ(d, r) *
        quadgk(x -> x^2(β-1) * ρ(d, x) * M(d, x), r, utilities.get_inf(r))[1])
end
# function σr(d::AbstractDensity, r::T, β::T) where T <: Real
#     value = sqrt(quadgk(x -> x^2β * ρ(d, x) * ∇Φ(d, x), r, utilities.get_inf(r))[1] / r^(2β) / ρ(d, r))
#     # Convert from km^(1/2) kpc^(1/2) / s to km/s
#     return value * 2947102009410051//16777216 # 1.7566096838772e8
# end
"""
    σlos(d::AbstractDensity, r::Real, β::Real)

Returns the line-of-sight projected velocity dispersion [km/s] of `d` at projected radius `r` [kpc] for constant velocity anisotropy `β` given by

```math
\\sigma_{\\text{LOS}}^2 \\left( R \\right) = \\frac{2}{\\Sigma \\left( R \\right)} \\int_R^\\infty \\left(1 - \\beta \\frac{R^2}{r^2} \\right) \\frac{r \\, \\rho \\left( r \\right) \\, \\sigma_r^2 \\left( r \\right)}{\\sqrt{r^2 - R^2}} \\, dr
```
"""
σlos(d::AbstractDensity, r::Real, β::Real) =
    sqrt(2 / Σ(d,r) * quadgk(x->(1-β*(r/x)^2) * ρ(d,x) * σr(d,x,β)^2 * x / sqrt(x^2 - r^2), r, utilities.get_inf(r))[1])
"""
    Φ(d::AbstractDensity, r::Real)
    Φ(uu::GalaxyProfiles.ΦdimensionUnits, d::AbstractDensity, r::Real)
    Φ(d::AbstractDensity, r::Unitful.Length)
    Φ(uu::GalaxyProfiles.ΦdimensionUnits, d::AbstractDensity, r::Unitful.Length)

Returns the potential [km^2 s^-2] of the density distribution `d` at radius `r` [kpc].

This is typically defined as
```math
\\begin{aligned}
\\Phi(R) &= -4\\pi G \\left( \\frac{1}{R} \\int_0^R r^2 \\rho(r) dr + \\int_R^\\infty r \\rho(r) dr \\right) \\newline
&= -G \\left( \\frac{M (R)}{R} - \\int_R^\\infty r \\rho(r) dr \\right) \\newline
% \\Phi(R) &= -\\frac{G}{R} \\int_0^R dM(r) - G \\int_R^\\infty \\frac{dM(r)}{r} \\newline
\\end{aligned}
```

where ``M(R)`` is the mass internal to radius ``R``. The final expression is used in the generic implementation of `Φ`. Given that there is also a generic implementation of [``M(R)``](@ref GalaxyProfiles.M) that calculates it from the density [`ρ`](@ref), the potential `Φ` can be calculated from a density profile by defining only the density [`ρ`](@ref).

If an analytic expression for ``M(r)`` exists, then an equivalent equation expressed in terms of the enclosed mass may be used,

```math
\\begin{aligned}
\\Phi(R) &= -4\\pi G \\left( \\frac{1}{R} \\int_0^R r^2 \\rho(r) dr + \\int_R^\\infty r \\rho(r) dr \\right) \\newline
&= -G \\left( \\frac{M (R)}{R} - \\int_R^\\infty \\frac{dM(r)}{r} \\right) \\newline
&= -G \\, \\int_R^\\infty \\frac{M \\left( r \\right)}{r^2} dr \\newline
\\end{aligned}
```

The above integrals are not finite for some mass distributions (e.g., [`GeneralIsothermal`](@ref) with some choices of `α`). In these cases, it is convention to define the potential ar `r` as the potential difference between `r` and the characteristic scale radius of the distribution; i.e.
```math
\\Phi(R) - \\Phi(R_s) = G \\int_{R_s}^R \\frac{M(r)}{r^2} dr
```
such that ``\\Phi(R_s)\\equiv0``.
"""
Φ(d::AbstractDensity{T}, r::S) where {T <: Real, S <: Real} =
    -convert(promote_type(T,S), constants.Gvelkpc) *
    (M(d,r)/r + quadgk(x->x * ρ(d,x), r, utilities.get_inf(T))[1] * 4 * π)
"""
    ∇Φ(d::AbstractDensity, r::Real)
    ∇Φ(uu::u.AccelerationUnits, d::AbstractDensity, r::Real)
    ∇Φ(d::AbstractDensity, r::Unitful.Length)
    ∇Φ(uu::u.AccelerationUnits, d::AbstractDensity, r::Unitful.Length)

Returns the radial derivative of the potential [`Φ(d,r)`](@ref Φ) [km s^-2] of profile `d` at radius `r` [kpc].

For a spherical mass distribution, the gravitational acceleration term in Newton's second law ``F=m*a`` is ``a = -G M(r) / r^2 = -\\nabla\\Phi`` by definition. This is derived from the definition of Φ below.

# Derivation from Φ

Recall that the potential can be written as

```math
\\Phi(R) = -G \\, \\int_R^\\infty \\frac{M \\left( r \\right)}{r^2} dr.
```

so that our derivative is

```math
\\begin{aligned}
\\nabla\\Phi(R) &= \\frac{\\partial \\Phi (R)}{\\partial R} \\newline
&= \\frac{\\partial}{\\partial R} \\left[ -G \\, \\int_R^\\infty \\frac{M \\left( r \\right)}{r^2} dr \\right] \\newline
\\end{aligned}
```

We are thus taking a radial derivative of a radial integral. By applying the fundamental thereom of calculus and remembering that we desire ``\\lim_{R \\to \\infty} \\Phi(R) = 0``, we can simply write

```math
\\begin{aligned}
\\nabla\\Phi(R) &= \\frac{G M \\left( r \\right)}{r^2}
\\end{aligned}
```
"""
∇Φ(d::AbstractDensity{T}, r::S) where {T <: Real, S <: Real} =
    convert(promote_type(T,S), constants.Gvelkpc2) * M(d, r) / r^2
"""
    ∇∇Φ(d::AbstractDensity, r::Real)
    ∇∇Φ(uu::GalaxyProfiles.∇∇ΦdimensionUnits, d::AbstractDensity, r::Real)
    ∇∇Φ(d::AbstractDensity, r::Unitful.Length)
    ∇∇Φ(uu::GalaxyProfiles.∇∇ΦdimensionUnits, d::AbstractDensity, r::Unitful.Length)

Returns the second radial derivative of the potential `Φ(d,r)` [km s^-2 kpc^-1] evaluated at radius `r` [kpc].

As the first radial derivative ``\\nabla\\Phi(R) = \\frac{G M \\left( r \\right)}{r^2}``, by the product rule the second radial derivative is

```math
\\nabla\\nabla\\Phi(R) = G \\left( \\frac{\\nabla M(R)}{r^2} - \\frac{2 M(R)}{r^3} \\right)
```
which is used as the generic implementation through the methods [`M`](@ref) and [`∇M`](@ref).

Note that this is *not* the same as the Laplacian operator that appears the Poisson equation ``\\nabla^2 \\Phi(R) = 4 \\pi G \\rho(R)``. In spherical coordinates, the radial component of the left hand side of Poisson's equation expands to

```math
\\nabla^2 \\Phi(R) = \\frac{1}{R^2} \\frac{\\partial}{\\partial R} \\left( R^2 \\frac{\\partial \\Phi(R)}{\\partial R} \\right)
```

which is not equivalent to the second radial gradient ``\\frac{\\partial^2 \\Phi(R)}{\\partial R^2}``, which is what this method returns.
"""
∇∇Φ(d::AbstractDensity{T}, r::S) where {T <: Real, S <: Real} =
    convert(promote_type(T, S), constants.Gvelkpc2) * (∇M(d, r) / r^2 - 2M(d, r) / r^3)

# check_vals(d::AbstractMassProfile,a::Number,T::DataType,S::DataType) = (promote_type(T,S), promote(params(d)...,a))
