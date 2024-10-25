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

Evaluate the density of `d` at radius `r`.
"""
ρ(d::AbstractDensity, r::Real)
"""
    invρ(d::AbstractMassProfile, x::Real
    invρ(d::AbstractMassProfile, x::Real,
        interval::NTuple{2,Real}=(scale_radius(d)/100,100*scale_radius(d)); kws...)
    invρ(uu::Unitful.LengthUnits, d::AbstractMassProfile, x::Real)
    invρ(d::AbstractMassProfile, x::Unitful.Density)
    invρ(uu::Unitful.LengthUnits, d::AbstractMassProfile, x::Unitful.Density)

Solve for the radius `r` at which the density is `x` for profile `d`. Requires `x>0`. When this method is not specialized for `d`, it will use an interval bracketing method from [`Roots.jl`](https://github.com/JuliaMath/Roots.jl), requiring that `ρ(d,r)` be defined. For this method, `kws...` are passed to `Roots.find_zero`.
"""
invρ(d::AbstractDensity, x::T, interval::NTuple{2,S}=(scale_radius(d)/100,100*scale_radius(d)); kws...) where {T<:Real, S<:Real} = (U = promote_type(T, S); x <= 0 ? throw(DomainError(x,"x must be > 0")) : find_zero(y->ρ(d,y)-x,(U(interval[1]),U(interval[2])); kws...))
"""
    ∇ρ(d::AbstractDensity, r::Real)
    ∇ρ(uu::GalaxyProfiles.∇ρdimensionUnits, d::AbstractDensity, r::Real)
    ∇ρ(d::AbstractDensity, r::Unitful.Length)
    ∇ρ(uu::GalaxyProfiles.∇ρdimensionUnits, d::AbstractDensity, r::Unitful.Length)

The gradient of `ρ(d,r)` evaluated at radius `r`.
"""
∇ρ(d::AbstractDensity, r::Real)
"""
    ρmean(d::AbstractDensity, r::Real)
    ρmean(uu::Unitful.DensityUnits, d::AbstractDensity, r::Real)
    ρmean(d::AbstractDensity, r::Unitful.Length)
    ρmean(uu::Unitful.DensityUnits, d::AbstractDensity, r::Unitful.Length)

The average density inside `r`; defaults to enclosed mass divided by volume; `M(d,r) * 3 / (4π*r^3)`. 
"""
ρmean(d::AbstractDensity, r::Real) = M(d,r) * 3 / 4 / π / r^3
"""
    invρmean(d::AbstractDensity, x::Real)
    invρmean(d::AbstractDensity, x::Real,
        interval::NTuple{2,Real}=(scale_radius(d)/100,100*scale_radius(d)); kws...)
    invρmean(uu::Unitful.LengthUnits, d::AbstractDensity, x::Real)
    invρmean(d::AbstractDensity, x::Unitful.Density)
    invρmean(uu::Unitful.LengthUnits, d::AbstractDensity, x::Unitful.Density)

Solve for the radius `r` inside which the average density is `x`. Requires `x>0`. When this method is not specialized for `d`, it will use an interval bracketing method from [`Roots.jl`](https://github.com/JuliaMath/Roots.jl), requiring that `ρmean(d,r)` or `M(d,r)` be defined. For this method, `kws...` are passed to `Roots.find_zero`.
"""
invρmean(d::AbstractDensity, x::T, interval::NTuple{2,S}=(scale_radius(d)/100,100*scale_radius(d)); kws...) where {T<:Real, S<:Real} = (U = promote_type(T, S); x <= 0 ? throw(DomainError(x,"x must be > 0")) : find_zero(y->ρmean(d,y)-x,(U(interval[1]), U(interval[2])); kws...))
"""
    Σ(d::AbstractMassProfile, r::Real)
    Σ(uu::GalaxyProfiles.SurfaceDensityUnits, d::AbstractMassProfile, r::Real)
    Σ(d::AbstractMassProfile, r::Unitful.Length)
    Σ(uu::GalaxyProfiles.SurfaceDensityUnits, d::AbstractMassProfile, r::Unitful.Length)

Evaluate the surface density of `d` profile at radius `r`. For 3D density profiles (i.e., [`AbstractDensity`](@ref GalaxyProfiles.AbstractDensity)), this will be the projected surface density, which (for spherical systems) is defined via the Abel integral 

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

The gradient of `Σ(d,r)` evaluated at radius `r`.
"""
∇Σ(d::AbstractMassProfile, r::Real)
"""
    Σmean(d::AbstractMassProfile, r::Real)
    Σmean(uu::GalaxyProfiles.SurfaceDensityUnits, d::AbstractMassProfile, r::Real)
    Σmean(d::AbstractMassProfile, r::Unitful.Length)
    Σmean(uu::GalaxyProfiles.SurfaceDensityUnits, d::AbstractMassProfile, r::Unitful.Length)

Evaluates the mean projected surface density inside the radius `r`; defaults to `Mproj(d::AbstractMassProfile,r::Real) / (π * r^2)`.
"""
Σmean(d::AbstractMassProfile, r::T) where T <: Real = Mproj(d,r) / r^2 / π
"""
    invΣ(d::AbstractMassProfile, x::Real)
    invΣ(d::AbstractMassProfile, x::Real,
        interval::NTuple{2,Real}=(scale_radius(d)/100,100*scale_radius(d)); kws...)
    invΣ(uu::Unitful.LengthUnits, d::AbstractMassProfile, x::Real)
    invΣ(d::AbstractMassProfile, r::Unitful.Length)
    invΣ(uu::Unitful.LengthUnits, d::AbstractMassProfile, r::Unitful.Length)

Solve for the radius `r` at which the surface density is `x` for profile `d`. Requires `x>0`. When this method is not specialized for `d`, it will use an interval bracketing method from [`Roots.jl`](https://github.com/JuliaMath/Roots.jl), requiring that `Σ(d,r)` be defined. For this method, `kws...` are passed to `Roots.find_zero`.
"""
invΣ(d::AbstractMassProfile, x::T, interval::NTuple{2,S}=(scale_radius(d)/100,100*scale_radius(d)); kws...) where {T<:Real, S<:Real} = (U = promote_type(T, S); x <= 0 ? throw(DomainError(x,"x must be > 0")) : find_zero(y->Σ(d,y)-x,(U(interval[1]), U(interval[2])); kws...))
"""
    M(d::AbstractDensity, r::Real; kws...)
    M(uu::Unitful.MassUnits, d::AbstractDensity, r::Real)
    M(d::AbstractDensity, r::Unitful.Length)
    M(uu::Unitful.MassUnits, d::AbstractDensity, r::Unitful.Length)

Evaluate the total mass enclosed within a radius `r` for the profile `d`. For spherical systems this is given by the integral

```math
M(\\lt R) = 4\\pi \\int_0^R r^2 ρ(r) dr
```

When this method is not specialized for `d`, it will compute the numerical integral using `QuadGK.quadgk`;  the provided `kws...` will be passed through. 
"""
M(d::AbstractDensity, r::Real; kws...) = quadgk(x->x^2 * ρ(d,x), zero(r), r)[1] * 4 * π
"""
    ∇M(d::AbstractDensity, r::Real)
    ∇M(uu::GalaxyProfiles.∇mdimensionUnits, d::AbstractDensity, r::Real)
    ∇M(d::AbstractDensity, r::Unitful.Length)
    ∇M(uu::GalaxyProfiles.∇mdimensionUnits, d::AbstractDensity, r::Unitful.Length)

The gradient of `M(d,r)` evaluated at radius `r`.
"""
∇M(d::AbstractMassProfile, r::Real)
"""
    invM(d::AbstractDensity, x::Real)
    invM(d::AbstractDensity, x::Real,
        interval::NTuple{2,Real}=(scale_radius(d)/100,100*scale_radius(d)); kws...)
    invM(uu::Unitful.LengthUnits, d::AbstractDensity, x::Real)
    invM(d::AbstractDensity, x::Unitful.Mass)
    invM(uu::Unitful.LengthUnits, d::AbstractDensity, x::Unitful.Mass)

Solve for the radius `r` at which the enclosed mass is `x` for profile `d`. Requires `x>0`. When this method is not specialized for `d`, it will use an interval bracketing method from [`Roots.jl`](https://github.com/JuliaMath/Roots.jl), requiring that `M(d,r)` be defined. For this method, `kws...` are passed to `Roots.find_zero`.
"""
invM(d::AbstractDensity, x::T, interval::NTuple{2,S}=(scale_radius(d)/100,100*scale_radius(d)); kws...) where {T<:Real, S<:Real} = (U = promote_type(T, S); x <= 0 ? throw(DomainError(x,"x must be > 0")) : find_zero(y->M(d,y)-x,(U(interval[1]), U(interval[2])); kws...))
"""
    Mtot(d::AbstractMassProfile)
    Mtot(uu::Unitful.MassUnits, d::AbstractMassProfile)

Evaluate the total mass of `d`. Should give same answer as `M(d,Inf)`. 
"""
Mtot(d::AbstractMassProfile)
"""
    Mproj(d::AbstractMassProfile, r::Real)
    Mproj(uu::Unitful.MassUnits, d::AbstractMassProfile, r::Real)
    Mproj(d::AbstractMassProfile, r::Unitful.Length)
    Mproj(uu::Unitful.MassUnits, d::AbstractMassProfile, r::Unitful.Length)

Evaluate the total line-of-sight projected mass enclosed within a radius `r` for the profile `d`. For spherical systems this is given by the integral
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

The gradient of `Mproj(d,r)` evaluated at radius `r`.
"""
∇Mproj(d::AbstractMassProfile, r::Real)
"""
    invMproj(d::AbstractMassProfile, x::Real)
    invMproj(d::AbstractMassProfile, x::Real,
        interval::NTuple{2,Real}=(scale_radius(d)/100,100*scale_radius(d)); kws...)
    invMproj(d::AbstractMassProfile, x::Unitful.Mass)
    invMproj(uu::Unitful.LengthUnits, d::AbstractMassProfile, x::Unitful.Mass)

Solve for the radius `r` at which the line-of-sight projected enclosed mass is `x` for profile `d`. Requires `x>0`. When this method is not specialized for `d`, it will use an interval bracketing method from [`Roots.jl`](https://github.com/JuliaMath/Roots.jl), requiring that `M(d,r)` be defined. For this method, `kws...` are passed to `Roots.find_zero`.
"""
invMproj(d::AbstractMassProfile, x::T, interval::NTuple{2,S}=(scale_radius(d)/100,100*scale_radius(d)); kws...) where {T<:Real, S<:Real} = (U = promote_type(T, S); x <= 0 ? throw(DomainError(x,"x must be > 0")) : find_zero(y->Mproj(d,y)-x,(U(interval[1]), U(interval[2])); kws...))
"""
    dynamical_time(d::Union{AbstractDensity, AbstractMassProfile}, r::Real)
    dynamical_time(d::Union{AbstractDensity, AbstractMassProfile}, r::Unitful.Length)
    dynamical_time(uu::Unitful.TimeUnits, d::Union{AbstractDensity, AbstractMassProfile}, r::Unitful.Length)

Returns the dynamical time at radius `r` in the provided density profile. The dynamical time is a characteristic timescale associated with orbits in potentials. This implementation specifically calculates Equation 2.39 in Binney & Tremaine Galactic Dynamics 2E, with the replacement of the standard density [`ρ`](@ref) for the average density interior to `r`, [`ρmean`](@ref), to better account for inhomogenous systems. 

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

where ``M \\left( r \\right)}`` is the mass enclosed inside radius ``r``, provided by the method [`M`](@ref). This equation is used when the argument `d::AbstractMassProfile`.

# Alternative Definitions
Note that some authors prefer to omit the factor of ``\\sqrt{ \\frac{1}{16} } = 1/4`` in the denominator of the first equation above; this is, for example, the convention used by galpy. Dynamical times thus defined will be four times larger than those calculated by this method. 
"""
function dynamical_time(d::AbstractDensity, r::T) where {T <: Real}
    # constants.Gkpc has time units of 1/s^2; 31557600 is seconds in year to return value in yr
    # return sqrt(inv(GalaxyProfiles.ρmean(d, r)) * 3 * π / 16 / T(constants.Gkpc)) / 31557600
    return sqrt(inv(GalaxyProfiles.ρmean(d, r)) * 3 * π) / (4 * 31557600 * sqrt(T(constants.Gkpc)))
end
dynamical_time(d::AbstractMassProfile, r::T) where {T <: Real} =
    π * sqrt(r^3 / M(d, r)) / (2 * 31557600 * sqrt(T(constants.Gkpc)))
    # π * sqrt( r^3 / 4 / T(constants.Gkpc) / M(d, r) ) / 31557600
"""
    cdf2D(d::AbstractMassProfile, r::Real)
    cdf2D(d::AbstractMassProfile, r::Unitful.Quantity)

Evaluate the cumulative distribution function of the profile `d` at `r` in two dimensions (i.e., along a line of sight). This is defined as `Σ(d,r)/Mtot(d)`. 
"""
cdf2D(d::AbstractMassProfile, r::Real) = Σ(d,r) / Mtot(d)
"""
    cdf3D(d::AbstractDensity, r::Real)
    cdf3D(d::AbstractDensity, r::Unitful.Quantity)

Evaluate the cumulative distribution function of the profile `d` at `r` in three dimensions. This is defined as `M(d,r)/Mtot(d)`. 
"""
cdf3D(d::AbstractDensity, r::Real) = M(d,r) / Mtot(d)
"""
    ccdf2D(d::AbstractMassProfile, r::Real)
    ccdf2D(d::AbstractMassProfile, r::Unitful.Quantity)

Evaluate the complementary cumulative distribution function of the profile `d` at `r` in two dimensions (i.e., along a line of sight). This is defined as `1 - cdf2D(d,r) = 1 - Σ(d,r)/Mtot(d)`. 
"""
ccdf2D(d::AbstractMassProfile, r::Real) = 1 - cdf2D(d,r)
"""
    ccdf3D(d::AbstractMassProfile, r::Real)
    ccdf3D(d::AbstractMassProfile, r::Unitful.Quantity)

Evaluate the complementary cumulative distribution function of the profile `d` at `r` in three dimensions. This is defined as `1 - cdf3D(d,r) = 1 - M(d,r)/Mtot(d)`. 
"""
ccdf3D(d::AbstractMassProfile, r::Real) = 1 - cdf3D(d,r)
"""
    quantile2D(d::AbstractMassProfile, r::Real)
    quantile2D(d::AbstractMassProfile, r::Unitful.Quantity)

Evaluate the inverse cumulative distribution function of the profile `d` at `r` in two dimensions (i.e., along a line of sight).
"""
quantile2D(d::AbstractMassProfile, r::Real)
"""
    quantile3D(d::AbstractDensity, r::Real)
    quantile3D(d::AbstractDensity, r::Unitful.Quantity)

Evaluate the inverse cumulative distribution function of the profile `d` at `r` in three dimensions.
"""
quantile3D(d::AbstractDensity, r::Real)
"""
    cquantile2D(d::AbstractMassProfile, r::Real)
    cquantile2D(d::AbstractMassProfile, r::Unitful.Quantity)

Evaluate the complementary quantile of the profile `d` at `r` in two dimensions (i.e., along a line of sight). This is defined as `quantile2D(d, 1-r)`.
"""
cquantile2D(d::AbstractMassProfile, r::Real)
"""
    cquantile3D(d::AbstractDensity, r::Real)
    cquantile3D(d::AbstractDensity, r::Unitful.Quantity)

Evaluate the complementary quantile of the profile `d` at `r` in three dimensions. This is defined as `quantile3D(d, 1-r)`.
"""
cquantile3D(d::AbstractDensity, r::Real)
"""
    Vcirc(d::AbstractDensity, r::Real)
    Vcirc(uu::Unitful.VelocityUnits, d::AbstractDensity, r::Real)
    Vcirc(d::AbstractDensity, r::Unitful.Length)
    Vcirc(uu::Unitful.VelocityUnits, d::AbstractDensity, r::Unitful.Length)

Evaluate the circular velocity at `r`, defined as the speed of a particle of insignificant mass in a circular orbit at radius `r`. This is calculated as

```math
v_c^2(r) = \\frac{G M(r)}{r} = r \\frac{d\\Phi}{dr} = r\\nabla\\Phi(r)
```

By default uses `G` in units such that if `rs` and `r` are in kpc, the velocity ends up in `km/s`. For example, for [`GeneralIsothermal`](@ref) we have `[G] = [kpc * km^2 / Msun / s^2]` so that the velocity ends up in `km/s`. Falls back to `sqrt( GalaxyProfiles.constants.Gvelkpc * M(d::AbstractDensity,r) / r)`.
"""
Vcirc(d::AbstractDensity, r::T) where T <: Real = sqrt(T(constants.Gvelkpc) * M(d, r) / r)
"""
    Vesc(d::AbstractDensity, r::Real)
    Vesc(uu::Unitful.VelocityUnits, d::AbstractDensity, r::Real)
    Vesc(d::AbstractDensity, r::Unitful.Length)
    Vesc(uu::Unitful.VelocityUnits, d::AbstractDensity, r::Unitful.Length)

Evaluate the escape velocity at `r`, calculated as
```math
v^2_{\\text{esc}}(r) = 2 |\\Phi(r)|
```
if ``\\Phi \\to 0`` for ``r \\to \\infty``; see the note for [`Φ`](@ref).

By default uses `G` in units such that if `rs` and `r` are in kpc, the velocity ends up in `km/s`. For example, for [`GeneralIsothermal`](@ref) we have `[G] = [kpc * km^2 / Msun / s^2]` so that the velocity ends up in `km/s`. 
"""
Vesc(d::AbstractDensity,r::Real) = sqrt( 2 * abs(Φ(d, r)) )
"""
    Vmax(d::AbstractDensity)
    Vmax(uu::Unitful.VelocityUnits,d::AbstractDensity)

Returns the maximum circular velocity of `d` in [km/s] and the corresponding radius in [kpc]. Can be found by solving

```math
    \\frac{d v_c(r)}{dr} = 0
```
for `r`, where ``v_c`` is the circular velocity, then evaluating the circular velocity at `r`. 
"""
Vmax(d::AbstractDensity,r::Real)
"""
    σr(d::AbstractDensity, r::Real, β::Real)
Returns the radial velocity dispersion in [km/s] of `d` at radius `r` for constant velocity anisotropy `β` given by Equation 4.216 in Binney & Tremaine Galactic Dynamics 2E,

```math
\\sigma_r^2 \\left( R \\right) = \\frac{1}{R^{2\\beta} \\, \\rho(R)} \\int_R^\\infty r^{2\\beta} \\, \\rho(r) \\, \\frac{d\\Phi}{dr} \\, dr
```

and as ``\\frac{d\\Phi}{dr} = -G M(r) / r^2`` we can alternatively write

```math
\\sigma_r^2 \\left( R \\right) = \\frac{G}{R^{2\\beta} \\, \\rho(R)} \\int_R^\\infty r^{2\\left( \\beta-1 \\right)} \\, \\rho(r) \\, M\\left( r \\right) \\, dr
```

which is the default fallback method as we expect ``M(r)`` to be more commonly available than ``\\frac{d\\Phi}{dr}``.
"""
σr(d::AbstractDensity, r::T, β::S) where {T <: Real, S <: Real} = σr(d, promote(r, β)...)
σr(d::AbstractDensity, r::T, β::T) where T <: Real =
    sqrt(T(constants.Gvelkpc) * quadgk(x -> x^2(β-1) * ρ(d, x) * M(d, x), r, utilities.get_inf(r))[1] / r^(2β) / ρ(d, r))
# function σr(d::AbstractDensity, r::T, β::T) where T <: Real
#     value = sqrt(quadgk(x -> x^2β * ρ(d, x) * ∇Φ(d, x), r, utilities.get_inf(r))[1] / r^(2β) / ρ(d, r))
#     # Convert from km^(1/2) kpc^(1/2) / s to km/s
#     return value * 2947102009410051//16777216 # 1.7566096838772e8
# end
"""
    σlos(d::AbstractDensity, r::Real, β::Real)

Returns the line-of-sight projected velocity dispersion in [km/s] of `d` at projected radius `r` for constant velocity anisotropy `β` given by

```math
\\sigma_{\\text{LOS}}^2 \\left( R \\right) = \\frac{2}{\\Sigma \\left( R \\right)} \\int_R^\\infty \\left(1 - \\beta \\frac{R^2}{r^2} \\right) \\frac{r \\, \\rho \\left( r \\right) \\, \\sigma_r^2 \\left( r \\right)}{\\sqrt{r^2 - R^2}} \\, dr
```
"""
σlos(d::AbstractDensity, r::T, β::S) where {T <: Real, S <: Real} = σlos(d, promote(r, β)...)
σlos(d::AbstractDensity, r::T, β::T) where T <: Real =
    sqrt(2 / Σ(d,r) * quadgk(x->(1-β*(r/x)^2) * ρ(d,x) * σr(d,x,β)^2 * x / sqrt(x^2 - r^2), r, utilities.get_inf(r))[1])
"""
    Φ(d::AbstractDensity, r::Real)
    Φ(uu::GalaxyProfiles.ΦdimensionUnits, d::AbstractDensity, r::Real)
    Φ(d::AbstractDensity, r::Unitful.Length)
    Φ(uu::GalaxyProfiles.ΦdimensionUnits, d::AbstractDensity, r::Unitful.Length)

Evaluate the potential of the density distribution at `r`. This is typically defined as
```math
\\Phi(R) = -\\frac{G}{R} \\int_0^R dM(r) - G \\int_R^\\infty \\frac{dM(r)}{r}
= -4\\pi G \\left( \\frac{1}{R} \\int_0^R r^2 \\rho(r) dr + \\int_R^\\infty r \\rho(r) dr \\right)
```
However, these integrals are not finite for some mass distributions (e.g., [`GeneralIsothermal`](@ref) with some choices of `α`); in these cases, it is convention to define the potential ar `r` as the potential difference between `r` and the characteristic scale radius of the distribution; i.e.
```math
\\Phi(R) - \\Phi(R_s) = G \\int_{R_s}^R \\frac{M(r)}{r^2} dr
```
such that ``\\Phi(R_s)\\equiv0``.
By default uses `G` in units such that if `rs` and `r` are in kpc, the potential ends up in `(km/s)^2`.
"""
Φ(d::AbstractDensity, r::T) where T <: Real =
    -T(constants.Gvelkpc) * (M(d,r)/r + quadgk(x->x * ρ(d,x), r, utilities.get_inf(T))[1] * 4 * π)
"""
    ∇Φ(d::AbstractDensity, r::Real)
    ∇Φ(uu::u.AccelerationUnits, d::AbstractDensity, r::Real)
    ∇Φ(d::AbstractDensity, r::Unitful.Length)
    ∇Φ(uu::u.AccelerationUnits, d::AbstractDensity, r::Unitful.Length)

The gradient of the potential `Φ(d,r)` evaluated at radius `r`.
By default uses `G` in units such that if the length units of `d.rs`, `r`, and `d.ρ0` are kpc, the derivative of the potential is returned as `[km/s^2]`.
"""
∇Φ(d::AbstractDensity,r::Real)
"""
    ∇∇Φ(d::AbstractDensity, r::Real)
    ∇∇Φ(uu::GalaxyProfiles.∇∇ΦdimensionUnits, d::AbstractDensity, r::Real)
    ∇∇Φ(d::AbstractDensity, r::Unitful.Length)
    ∇∇Φ(uu::GalaxyProfiles.∇∇ΦdimensionUnits, d::AbstractDensity, r::Unitful.Length)

The second order gradient of the potential `Φ(d,r)` evaluated at radius `r`.
By default uses `G` in units such that if the length units of `d.rs`, `r`, and `d.ρ0` are kpc, the derivative of the potential is returned as `[km/s^2/kpc]`.
"""
∇∇Φ(d::AbstractDensity,r::Real)


# check_vals(d::AbstractMassProfile,a::Number,T::DataType,S::DataType) = (promote_type(T,S), promote(params(d)...,a))
