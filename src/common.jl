module constants
# const Gkpc = 4.5171030498949647e-39     # G in kpc^3/solMass/s^2
# const Gpc = 4.5171030498949647e-30      # G in pc^3/solMass/s^2
# const Gvelkpc = 4.30091727003628e-6     # G in kpc*km^2/solMass/s^2, for velocity calculations
# const Gvelkpc2 = 1.3938323614347172e-22 # G in kpc^2*km/solMass/s^2, for velocity calculations
# const Gvelpc = 4.30091727003628e-3      # G in kpc*km^2/solMass/s^2
const Gkpc = big"4.517103049894965029489870560544420380775973240960280169677794680278054678900598e-39" # G in kpc^3/solMass/s^2
const Gpc = big"4.517103049894965632856118045698970338735168159241032233763499248446748879359808e-30"  # G in pc^3/solMass/s^2
const Gvelkpc = big"4.3009172700362805113763897679746150970458984375e-6"   # G in kpc*km^2/solMass/s^2, for velocity and Φ calculations; (G |> ua.kpc*u.km^2/ua.Msun/u.s^2)
const Gvelkpc2 = big"1.393832361434717359905559114367183944790625806098294248158708796836435794830327e-22" # G in kpc^2*km/solMass/s^2, for ∇Φ calculations; (G |> ua.kpc^2*u.km/ua.Msun/u.s^2)
const Gvelpc = big"4.3009172700362805113763897679746150970458984375e-6"
const nfwvmaxpar = big"2.162581587064609834856553669603264573507154023820515448359231774428038693093796" # constant for NFW's Vmax
end
"""
    params(d::AbstractMassProfile)

Get the fields of the struct as expected by its defined methods.
"""
params(d::AbstractMassProfile)
"""
    scale_radius(d::AbstractMassProfile)
    scale_radius(uu::Unitful.LengthUnits,d::AbstractMassProfile)

Returns the characteristic scale radius of the profile; used for some default methods. An example is `rs` for the [`ExponentialDisk`](@ref) model. You should generally use `params` instead.
"""
scale_radius(d::AbstractMassProfile)
"""
    ρ(d::AbstractDensity, r::Real)
    ρ(uu::Unitful.DensityUnits, d::AbstractDensity, r::Real)
    ρ(d::AbstractDensity, r::Unitful.Length)
    ρ(uu::Unitful.DensityUnits, d::AbstractDensity, r::Unitful.Length)

Evaluate the density of `d` at radius `r`.
"""
ρ(d::AbstractDensity,r::Real)
"""
    invρ(d::AbstractMassProfile, x::Real
    invρ(d::AbstractMassProfile, x::Real, interval::Tuple=(scale_radius(d)/100,100*scale_radius(d))
    invρ(uu::Unitful.LengthUnits, d::AbstractMassProfile, x::Real)
    invρ(d::AbstractMassProfile, x::Unitful.Density)
    invρ(uu::Unitful.LengthUnits, d::AbstractMassProfile, x::Unitful.Density)

Solve for the radius `r` at which the density is `x` for profile `d`. Requires `x>0`. When this method is not specialized for `d`, it will use an interval bracketing method from [`Roots.jl`](https://github.com/JuliaMath/Roots.jl), requiring that `ρ(d,r)` be defined.
"""
invρ(d::AbstractDensity,x::Real,interval::Tuple=(scale_radius(d)/100,100*scale_radius(d))) = x <= 0 ? throw(DomainError(x,"x must be > 0")) : find_zero(y->ρ(d,y)-x,interval)
"""
    ∇ρ(d::AbstractDensity, r::Real)
    ∇ρ(uu::GalaxyProfiles.∇ρdimensionUnits, d::AbstractDensity, r::Real)
    ∇ρ(d::AbstractDensity, r::Unitful.Length)
    ∇ρ(uu::GalaxyProfiles.∇ρdimensionUnits, d::AbstractDensity, r::Unitful.Length)

The gradient of `ρ(d,r)` evaluated at radius `r`.
"""
∇ρ(d::AbstractDensity,r::Real)
"""
    ρmean(d::AbstractDensity, r::Real)
    ρmean(uu::Unitful.DensityUnits, d::AbstractDensity,r::Real)
    ρmean(d::AbstractDensity, r::Unitful.Length)
    ρmean(uu::Unitful.DensityUnits, d::AbstractDensity, r::Unitful.Length)

The average density inside `r`; defaults to enclosed mass divided by volume; `M(d,r) * 3 / (4π*r^3)`. 
"""
ρmean(d::AbstractDensity,r::Real) = M(d,r) * 3 / (4π*r^3)
"""
    invρmean(d::AbstractDensity, x::Real)
    invρmean(d::AbstractDensity, x::Real, interval::Tuple=(scale_radius(d)/100,100*scale_radius(d)))
    invρmean(uu::Unitful.LengthUnits, d::AbstractDensity, x::Real)
    invρmean(d::AbstractDensity, x::Unitful.Density)
    invρmean(uu::Unitful.LengthUnits, d::AbstractDensity, x::Unitful.Density)

Solve for the radius `r` inside which the average density is `x`. Requires `x>0`. When this method is not specialized for `d`, it will use an interval bracketing method from [`Roots.jl`](https://github.com/JuliaMath/Roots.jl), requiring that `ρmean(d,r)` or `M(d,r)` be defined.
"""
invρmean(d::AbstractDensity,x::Real,interval::Tuple=(scale_radius(d)/100,100*scale_radius(d))) = x <= 0 ? throw(DomainError(x,"x must be > 0")) : find_zero(y->ρmean(d,y)-x,interval)
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
Σ(d::AbstractMassProfile,r::Real) = 2 * quadgk(x->x*ρ(d,x)/sqrt(x^2-r^2),r,Inf)[1]
"""
    ∇Σ(d::AbstractMassProfile, r::Real)
    ∇Σ(uu::Unitful.DensityUnits, d::AbstractMassProfile, r::Real)
    ∇Σ(d::AbstractMassProfile, r::Unitful.Length)
    ∇Σ(uu::Unitful.DensityUnits, d::AbstractMassProfile, r::Unitful.Length)

The gradient of `Σ(d,r)` evaluated at radius `r`.
"""
∇Σ(d::AbstractMassProfile,r::Real)
"""
    Σmean(d::AbstractMassProfile, r::Real)
    Σmean(uu::GalaxyProfiles.SurfaceDensityUnits, d::AbstractMassProfile, r::Real)
    Σmean(d::AbstractMassProfile, r::Unitful.Length)
    Σmean(uu::GalaxyProfiles.SurfaceDensityUnits, d::AbstractMassProfile, r::Unitful.Length)

Evaluates the mean projected surface density inside the radius `r`; defaults to `Mproj(d::AbstractMassProfile,r::Real) / (π * r^2)`.
"""
Σmean(d::AbstractMassProfile,r::Real) = Mproj(d,r) / (π * r^2)
"""
    invΣ(d::AbstractMassProfile, x::Real)
    invΣ(d::AbstractMassProfile, x::Real, interval::Tuple=(scale_radius(d)/100,100*scale_radius(d))
    invΣ(uu::Unitful.LengthUnits, d::AbstractMassProfile, x::Real)
    invΣ(d::AbstractMassProfile, r::Unitful.Length)
    invΣ(uu::Unitful.LengthUnits, d::AbstractMassProfile, r::Unitful.Length)

Solve for the radius `r` at which the surface density is `x` for profile `d`. Requires `x>0`. When this method is not specialized for `d`, it will use an interval bracketing method from [`Roots.jl`](https://github.com/JuliaMath/Roots.jl), requiring that `Σ(d,r)` be defined.
"""
invΣ(d::AbstractMassProfile,x::Real,interval::Tuple=(scale_radius(d)/100,100*scale_radius(d))) = x <= 0 ? throw(DomainError(x,"x must be > 0")) : find_zero(y->Σ(d,y)-x,interval)
"""
    M(d::AbstractMassProfile, r::Real)
    M(uu::Unitful.MassUnits, d::AbstractMassProfile, r::Real)
    M(d::AbstractMassProfile, r::Unitful.Length)
    M(uu::Unitful.MassUnits, d::AbstractMassProfile, r::Unitful.Length)

Evaluate the total mass enclosed within a radius `r` for the profile `d`. For spherical systems this is given by the integral
```math
M(\\lt R) = 4\\pi \\int_0^R r^2 ρ(r) dr
```
"""
M(d::AbstractMassProfile,r::Real)
"""
    ∇M(d::AbstractMassProfile, r::Real)
    ∇M(uu::GalaxyProfiles.∇mdimensionUnits, d::AbstractMassProfile, r::Real)
    ∇M(d::AbstractMassProfile, r::Unitful.Length)
    ∇M(uu::GalaxyProfiles.∇mdimensionUnits, d::AbstractMassProfile, r::Unitful.Length)

The gradient of `M(d,r)` evaluated at radius `r`.
"""
∇M(d::AbstractMassProfile,r::Real)
"""
    invM(d::AbstractMassProfile, x::Real)
    invM(d::AbstractMassProfile, x::Real, interval::Tuple=(scale_radius(d)/100,100*scale_radius(d))
    invM(uu::Unitful.LengthUnits, d::AbstractMassProfile, x::Real)
    invM(d::AbstractMassProfile, x::Unitful.Mass)
    invM(uu::Unitful.LengthUnits, d::AbstractMassProfile, x::Unitful.Mass)

Solve for the radius `r` at which the enclosed mass is `x` for profile `d`. Requires `x>0`. When this method is not specialized for `d`, it will use an interval bracketing method from [`Roots.jl`](https://github.com/JuliaMath/Roots.jl), requiring that `M(d,r)` be defined.
"""
invM(d::AbstractMassProfile,x::Real,interval::Tuple=(scale_radius(d)/100,100*scale_radius(d))) = x <= 0 ? throw(DomainError(x,"x must be > 0")) : find_zero(y->M(d,y)-x,interval)
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
Mproj(d::AbstractMassProfile,r::Real) = 2π * quadgk(x->x * Σ(d,x),0,r)[1]
"""
    ∇Mproj(d::AbstractMassProfile, r::Real)
    ∇Mproj(uu::GalaxyProfiles.∇mdimensionUnits, d::AbstractMassProfile, r::Real)
    ∇Mproj(d::AbstractMassProfile, r::Unitful.Length)
    ∇Mproj(uu::GalaxyProfiles.∇mdimensionUnits, d::AbstractMassProfile, r::Unitful.Length)

The gradient of `Mproj(d,r)` evaluated at radius `r`.
"""
∇Mproj(d::AbstractMassProfile,r::Real)
"""
    invMproj(d::AbstractMassProfile, x::Real)
    invMproj(d::AbstractMassProfile, x::Real, interval::Tuple=(scale_radius(d)/100,100*scale_radius(d))
    invMproj(d::AbstractMassProfile, x::Unitful.Mass)
    invMproj(uu::Unitful.LengthUnits, d::AbstractMassProfile, x::Unitful.Mass)

Solve for the radius `r` at which the line-of-sight projected enclosed mass is `x` for profile `d`. Requires `x>0`. When this method is not specialized for `d`, it will use an interval bracketing method from [`Roots.jl`](https://github.com/JuliaMath/Roots.jl), requiring that `M(d,r)` be defined.
"""
invMproj(d::AbstractMassProfile,x::Real,interval::Tuple=(scale_radius(d)/100,100*scale_radius(d))) = x <= 0 ? throw(DomainError(x,"x must be > 0")) : find_zero(y->Mproj(d,y)-x,interval)
"""
    cdf(d::AbstractMassProfile, r::Real)
    cdf(d::AbstractMassProfile, r::Unitful.Quantity)

Evaluate the cumulative distribution function of the profile `d` at `r`. This is defined as `M(d,r)/Mtot(d)`. 
"""
cdf(d::AbstractMassProfile,r::Real) = M(d,r) / Mtot(d)
"""
    ccdf(d::AbstractMassProfile, r::Real)
    ccdf(d::AbstractMassProfile, r::Unitful.Quantity)

Evaluate the complementary cumulative distribution function of the profile `d` at `r`. This is defined as `1 - cdf(d,r) = 1 - M(d,r)/Mtot(d)`. 
"""
ccdf(d::AbstractMassProfile,r::Real) = 1 - cdf(d,r)
"""
    quantile(d::AbstractMassProfile, r::Real)
    quantile(d::AbstractMassProfile, r::Unitful.Quantity)

Evaluate the inverse cumulative distribution function of the profile `d` at `r`. 
"""
quantile(d::AbstractMassProfile,r::Real)
"""
    cquantile(d::AbstractMassProfile, r::Real)
    cquantile(d::AbstractMassProfile, r::Unitful.Quantity)

Evaluate the complementary quantile (i.e., `quantile(d,1-r)`) of the profile `d` at `r`. 
"""
cquantile(d::AbstractMassProfile,r::Real)
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
Vcirc(d::AbstractDensity,r::T) where {T<:Real} = sqrt( T(constants.Gvelkpc) * M(d,r) / r)
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
Vesc(d::AbstractDensity,r::Real) = sqrt( 2 * abs(Φ(d,r)) )
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
Φ(d::AbstractDensity,r::Real)
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
