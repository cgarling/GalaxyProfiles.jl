
"""
    ExponentialDisk(rs::Real,ρ0::Real)
    ExponentialDisk(rs::Real;ρ0=nothing,M=nothing)

Type descibring isotropic exponential surface density profiles with central surface density `ρ0` and scale radius `rs`. The surface density profile is

```math
\\Sigma(r) = \\rho_0 \\times \\exp \\left( \\frac{-r}{R_s} \\right)
```
"""
struct ExponentialDisk{T<:Real} <: AbstractSurfaceDensity
    Σ0::T
    rs::T
    # ExponentialDisk{T}(ρ0::T,rs::T) where {T} = new{T}(μ,σ)
end
# ExponentialDisk(rs::T;M=nothing,rs=nothing)
# ExponentialDisk(ρ0::T,rs::T) = 
ExponentialDisk(rs::Real,Σ0::Real) = ExponentialDisk(promote(rs,Σ0)...)
exponential_disk_from_M(M::T,rs::T) where {T<:Real} = ExponentialDisk(M/(8π*rs^3), rs)
exponential_disk_from_M(M::Real,rs::Real) = exponential_disk_from_M(promote(M,rs)...)
ExponentialDisk(rs::Real;M=nothing,Σ0=nothing) = isnothing(M) ? (@assert !isnothing(Σ0); ExponentialDisk(Σ0,rs)) : exponential_disk_from_M(M,rs)

Σ(d::ExponentialDisk,r::Real) = d.Σ0 * exp(-r/d.rs)
