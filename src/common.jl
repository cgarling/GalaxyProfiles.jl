# module constants
# const twopi = 2π
# end
"""
    params(d::AbstractMassProfile)

Get the fields of the struct as expected by its defined methods.
"""
params(d::AbstractMassProfile)
"""
    scale_radius(d::AbstractMassProfile)

Returns the characteristic scale radius of the profile; used for some default methods. An example is `rs` for the [`ExponentialDisk`](@ref) model. You should generally use `params` instead.
"""
scale_radius(d::AbstractMassProfile)
"""
    ρ(d::AbstractDensity, r::Real)
    ρ(d::AbstractDensity, r::Unitful.Quantity)

Evaluate the density of `d` at radius `r`.
"""
ρ(d::AbstractDensity,r::Real)
"""
    dρ_dr(d::AbstractDensity, r::Real)
    dρ_dr(d::AbstractDensity, r::Unitful.Quantity)

Evaluate the derivative of `ρ(d,r)` with respect to radius at `r`.
"""
dρ_dr(d::AbstractDensity,r::Real)
"""
    invρ(d::AbstractMassProfile, x::Real)
    invρ(d::AbstractMassProfile, r::Unitful.Quantity)
    invρ(d::AbstractMassProfile, x::Real, interval::Tuple=(scale_radius(d)/100,100*scale_radius(d))

Solve for the radius `r` at which the density is `x` for profile `d`. Requires `x>0`. When this method is not specialized for `d`, it will use an interval bracketing method from [`Roots.jl`](https://github.com/JuliaMath/Roots.jl), requiring that `ρ(d,r)` be defined.
"""
invρ(d::AbstractMassProfile,x::Real,interval::Tuple=(scale_radius(d)/100,100*scale_radius(d))) = x <= 0 ? throw(DomainError(x,"x must be > 0")) : find_zero(y->ρ(d,y)-x,interval)
"""
    Σ(d::AbstractMassProfile, r::Real)
    Σ(d::AbstractMassProfile, r::Unitful.Quantity)

Evaluate the surface density of `d` profile at radius `r`. For 3D density profiles (i.e., [`AbstractDensity`](@ref GalaxyProfiles.AbstractDensity)), this will be the projected surface density, which (for spherical systems) is defined via the Abel integral 

```math
\\Sigma(r) = 2 \\int_R^\\infty \\rho(r) \\frac{r}{ \\sqrt{r^2-R^2} } dr
```

which has an inverse of

```math
\\rho(r) = -\\frac{1}{\\pi} \\int_r^\\infty \\frac{d\\Sigma(R)}{dR} \\frac{dR}{\\sqrt{R^2-r^2}}
```
"""
Σ(d::AbstractMassProfile,r::Real)
"""
    dΣ_dr(d::AbstractMassProfile, r::Real)
    dΣ_dr(d::AbstractMassProfile, r::Unitful.Quantity)

Evaluate the derivative of the surface density `Σ(d,r)` with respect to radius of `d` at radius `r`.
"""
dΣ_dr(d::AbstractMassProfile,r::Real)
"""
    invΣ(d::AbstractMassProfile, x::Real)
    invΣ(d::AbstractMassProfile, r::Unitful.Quantity)
    invΣ(d::AbstractMassProfile, x::Real, interval::Tuple=(scale_radius(d)/100,100*scale_radius(d))

Solve for the radius `r` at which the surface density is `x` for profile `d`. Requires `x>0`. When this method is not specialized for `d`, it will use an interval bracketing method from [`Roots.jl`](https://github.com/JuliaMath/Roots.jl), requiring that `Σ(d,r)` be defined.
"""
invΣ(d::AbstractMassProfile,x::Real,interval::Tuple=(scale_radius(d)/100,100*scale_radius(d))) = x <= 0 ? throw(DomainError(x,"x must be > 0")) : find_zero(y->Σ(d,y)-x,interval)
"""
    M(d::AbstractMassProfile, r::Real)
    M(d::AbstractMassProfile, r::Unitful.Quantity)

Evaluate the total mass enclosed within a radius `r` for the profile `d`.
"""
M(d::AbstractMassProfile,r::Real)
"""
    dM_dr(d::AbstractMassProfile, r::Real)
    dM_dr(d::AbstractMassProfile, r::Unitful.Quantity)

Evaluate the derivative of the mass enclosed within a radius `r` with respect to `r` for the profile `d`.
"""
dM_dr(d::AbstractMassProfile,r::Real)
"""
    invM(d::AbstractMassProfile, x::Real)
    invM(d::AbstractMassProfile, r::Unitful.Quantity)
    invM(d::AbstractMassProfile, x::Real, interval::Tuple=(scale_radius(d)/100,100*scale_radius(d))

Solve for the radius `r` at which the enclosed mass is `x` for profile `d`. Requires `x>0`. When this method is not specialized for `d`, it will use an interval bracketing method from [`Roots.jl`](https://github.com/JuliaMath/Roots.jl), requiring that `M(d,r)` be defined.
"""
invM(d::AbstractMassProfile,x::Real,interval::Tuple=(scale_radius(d)/100,100*scale_radius(d))) = x <= 0 ? throw(DomainError(x,"x must be > 0")) : find_zero(y->M(d,y)-x,interval)
"""
    Mtot(d::AbstractMassProfile, r::Real)
    Mtot(d::AbstractMassProfile, r::Unitful.Quantity)

Evaluate the total mass of `d`. Should give same answer as `M(d,Inf)`. 
"""
Mtot(d::AbstractMassProfile,r::Real)
"""
    cdf(d::AbstractMassProfile, r::Real)
    cdf(d::AbstractMassProfile, r::Unitful.Quantity)

Evaluate the cumulative distribution function of the profile `d` at `r`. This is defined as `M(d,r)/Mtot(d)`. 
"""
cdf(d::AbstractMassProfile,r::Real)
"""
    ccdf(d::AbstractMassProfile, r::Real)
    ccdf(d::AbstractMassProfile, r::Unitful.Quantity)

Evaluate the complementary cumulative distribution function of the profile `d` at `r`. This is defined as `1 - cdf(d,r) = 1 - M(d,r)/Mtot(d)`. 
"""
ccdf(d::AbstractMassProfile,r::Real)
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

