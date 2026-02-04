
"""
    ExponentialDisk(Σ0::Real,rs::Real)
    ExponentialDisk(rs::Real;Σ0=nothing,M=nothing)

Type describing projected isotropic exponential surface density profiles with central surface density `Σ0` and scale radius `rs`. The surface density profile is

```math
\\Sigma(r) = \\rho_0 \\times \\exp \\left( \\frac{-r}{R_s} \\right)
```

The fields of `ExponentialDisk` are `Σ0, rs`. There are no methods defined for `ExponentialDisk` which use physical constants with units (e.g., `G`), so as long as `ExponentialDisk.rs` and the radius `r` you provide to methods are in the same units, and `Σ0` is in units of `[M/[r,rs]^2]`everything will work out. Generally just want to make sure the length units are uniform.

The following public methods are defined on this type:
 - [`Σ`](@ref), [`∇Σ`](@ref), [`invΣ`](@ref), [`Mproj`](@ref), [`∇Mproj`](@ref), [`invMproj`](@ref), [`Mtot`](@ref), [`cdf2D`](@ref), [`ccdf2D`](@ref), [`quantile2D`](@ref), [`cquantile2D`](@ref)

# See also
 - Convenience constructor [`ExponentialDiskDHI`](@ref).
"""
struct ExponentialDisk{TΣ, Tr} <: AbstractSurfaceDensity
    Σ0::TΣ
    rs::Tr
end
exponential_disk_from_M(M, rs) = ExponentialDisk(M/(2*π*rs^2), rs)
ExponentialDisk(rs; M=nothing, Σ0=nothing) = isnothing(M) ? (@assert !isnothing(Σ0); ExponentialDisk(Σ0,rs)) : exponential_disk_from_M(M, rs)
"""
    ExponentialDiskDHI(DHI::Real, MHI::Real, ΣDHI::Real=10^6)
    ExponentialDiskDHI(DHI::Unitful.Length, MHI::Unitful.Mass,
        ΣDHI::GalaxyProfiles.SurfaceDensity=1*UnitfulAstro.Msun/UnitfulAstro.pc^2)

Convenience constructor for the [`ExponentialDisk`](@ref) surface density profile. This will take a diameter for the disk `DHI` (e.g., the [diameter of a neutral hydrogen disk](https://ui.adsabs.harvard.edu/abs/2016MNRAS.460.2143W/abstract)), its total mass `MHI`, and `ΣDHI`, the surface density at `DHI`, and return an [`ExponentialDisk`](@ref) object with the correct central density and scale radius.

By default, `ΣDHI` is 10^6 solar masses per square kiloparsec (equivalent to one solar mass per square parsec), such that you should provide `DHI` in kpc if you are using `Real` inputs.

## Notes
At high masses (e.g., ``10^{10} \\ \\text{M}_\\odot < \\text{M}_{\\text{HI}}``) the numerical inversion is impossible; in this case, the exponential disk scale radius is set to `DHI/4`. 
"""
function ExponentialDiskDHI(DHI, MHI, ΣDHI=10^6) # DHI in kpc, M in Msun, ΣDHI in Msun/kpc^2
    # C = -0.5 * DHI * sqrt(π * ΣDHI / (2*MHI))
    C = -1//2 * DHI * sqrt(π * ΣDHI / 2 / MHI)
    if C > -1/ℯ
        rs = -DHI/(4 * lambertw(C,0))
    else
        rs = DHI/4
    end
    Σ0 = MHI / 2 / π / rs^2 # Σ0 = MHI/(2π * rs^2)
    return ExponentialDisk(Σ0,rs)
end

#### Parameters
"""
    (d.Σ0, d.rs) = params(d::ExponentialDisk)
"""
params(d::ExponentialDisk) = (d.Σ0,d.rs)
"""
    d.rs = scale_radius(d::ExponentialDisk)
"""
scale_radius(d::ExponentialDisk) = d.rs

#### Evaluation
function Σ(d::ExponentialDisk, r::Real) 
    Σ0, rs = params(d)
    return Σ0 * exp(-r/rs)
end
function ∇Σ(d::ExponentialDisk, r::Real)
    Σ0, rs = params(d)
    return -Σ0 * exp(-r/rs) / rs
end
# Σmean
function invΣ(d::ExponentialDisk, x::Real)
    Σ0, rs = params(d)
    x<=0 ? throw(DomainError(x, "x must be greater 0")) : rs * log(Σ0/x)
end
function Mproj(d::ExponentialDisk, r::Real)
    Σ0, rs = params(d)
    return isinf(r) ? Mtot(d) : Σ0 * rs^2 * (1 - exp(-r/rs)*(r/rs + 1)) * 2 * π
end
function ∇Mproj(d::ExponentialDisk, r::Real)
    Σ0, rs = params(d)
    return r * Σ0 * exp(-r/rs) * 2 * π
end
function invMproj(d::ExponentialDisk, x::Real)
    Σ0, rs = params(d)
    return -rs * (1 + lambertw( ((x/Mtot(d)) - 1) / ℯ , -1) )
end
function Mtot(d::ExponentialDisk)
    Σ0, rs = params(d)
    return Σ0 * rs^2 * 2 * π
end
function ccdf2D(d::ExponentialDisk, r::Real)
    Σ0, rs = params(d)
    ee = exp(-r/rs)
    return r*ee/rs + ee
end
cdf2D(d::ExponentialDisk, r::Real) = 1 - ccdf2D(d, r)
function cquantile2D(d::ExponentialDisk, x::Real)
    Σ0, rs = params(d)
    return -rs * (1+lambertw(-x/ℯ,-1))
end
quantile2D(d::ExponentialDisk, x::Real) = cquantile2D(d, 1-x)
