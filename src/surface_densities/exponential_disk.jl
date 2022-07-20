
"""
    ExponentialDisk(Σ0::Real,rs::Real)
    ExponentialDisk(rs::Real;Σ0=nothing,M=nothing)

Type descibring isotropic exponential surface density profiles with central surface density `ρ0` and scale radius `rs`. The surface density profile is

```math
\\Sigma(r) = \\rho_0 \\times \\exp \\left( \\frac{-r}{R_s} \\right)
```

The fields of `ExponentialDisk` are `Σ0,rs`.
"""
struct ExponentialDisk{T<:Real} <: AbstractSurfaceDensity
    Σ0::T
    rs::T
    # ExponentialDisk{T}(ρ0::T,rs::T) where {T} = new{T}(μ,σ)
end
ExponentialDisk(Σ0::Real,rs::Real) = ExponentialDisk(promote(Σ0,rs)...)
exponential_disk_from_M(M::T,rs::T) where {T<:Real} = ExponentialDisk(M/(8π*rs^3), rs)
exponential_disk_from_M(M::Real,rs::Real) = exponential_disk_from_M(promote(M,rs)...)
ExponentialDisk(rs::Real;M=nothing,Σ0=nothing) = isnothing(M) ? (@assert !isnothing(Σ0); ExponentialDisk(Σ0,rs)) : exponential_disk_from_M(M,rs)

#### Parameters
"""
    (d.Σ0,d.rs) = params(d::ExponentialDisk)
"""
params(d::ExponentialDisk) = (d.Σ0,d.rs)
"""
    d.rs = scale_radius(d::ExponentialDisk)
"""
scale_radius(d::ExponentialDisk) = d.rs

#### Evaluation

Σ(d::ExponentialDisk,r::Real) = d.Σ0 * exp(-r/d.rs)
function invΣ(d::ExponentialDisk,x::Real)
    Σ0,rs = params(d)
    x<=0 ? throw(DomainError(x,"x must be greater 0")) : rs * log(Σ0/x)
end
function dΣ_dr(d::ExponentialDisk,r::Real)
    Σ0,rs=params(d)
    2π * r * Σ0 * exp(-r/rs) 
end
function M(d::ExponentialDisk,r::Real)
    Σ0,rs=params(d)
    ee = exp(-r/rs)
    # 2π * rs * Σ0 * (rs - r*exp(-r/rs) - rs*exp(-r/rs))
    2π * rs * Σ0 * (rs - r*ee - rs*ee)
end
function dM_dr(d::ExponentialDisk,r::Real)
    Σ0,rs=params(d)
    2π * r * Σ0 * exp(-r/rs) 
end
function Mtot(d::ExponentialDisk)
    Σ0,rs=params(d)
    2π * Σ0 * rs^2
end
function ccdf(d::ExponentialDisk,r::Real)
    Σ0,rs=params(d)
    ee = exp(-r/rs)
    r*ee/rs + ee
end
cdf(d::ExponentialDisk,r::Real) = 1-ccdf(d,r)
function cquantile(d::ExponentialDisk,x::Real)
    Σ0,rs=params(d)
    -rs * (1+lambertw(-x/ℯ,-1))
end
quantile(d::ExponentialDisk,x::Real) = cquantile(d,1-x)
