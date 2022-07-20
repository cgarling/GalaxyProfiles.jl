
"""
    GeneralIsothermal(ρ0::Real,rs::Real,α::Real)
    GeneralIsothermal(rs::Real,α::Real,M::Real,Rmax::Real)

Type describing general isothermal density profiles with scale radius `rs`, power-law index `α`, and central density `ρ0`. The surface density profile is

```math
\\rho(r) = \\rho_0 \\times \\left( \\frac{r}{R_s} \\right)^{-\\alpha}
```

The fields of `GeneralIsothermal` are `ρ0,rs,α`.

# See also
 - [`SIS`](@ref)
"""
struct GeneralIsothermal{T<:Real} <: AbstractSurfaceDensity
    ρ0::T
    rs::T
    α::T
    # ExponentialDisk{T}(ρ0::T,rs::T) where {T} = new{T}(μ,σ)
end
GeneralIsothermal(ρ0::Real,rs::Real,α::Real) = ExponentialDisk(promote(ρ0,rs,α)...)
generalisothermal_from_M_Rmax(rs::T,α::T,M::T,Rmax::T) where {T<:Real} = α > 3 ? throw(DomainError(α,"Enclosed mass is only finite for α<3.")) : GeneralIsothermal((3 - α) * M / (4π * rs^α * Rmax^(3-α)),rs,α)
# generalisothermal_from_M_Rmax(rs::Real,α::Real,M::Real,Rmax::Real) = generalisothermal_from_M_Rmax(promote(rs,α,M,Rmax)...)
GeneralIsothermal(rs::Real,α::Real,M::Real,Rmax::Real) = generalisothermal_from_M_Rmax(promote(rs,α,M,Rmax)...)

"""
    SIS(ρ0::Real,rs::Real)
    SIS(rs::Real,M::Real,Rmax::Real)

Convenience function to construct a singular isothermal sphere; e.g. a [`GeneralIsothermal`](@ref) with `α=2`.
"""
SIS(ρ0::Real,rs::Real) = GeneralIsothermal(ρ0,rs,2)
SIS(rs::Real,M::Real,Rmax::Real) = GeneralIsothermal(rs,2,M,Rmax)

#### Parameters

params(d::GeneralIsothermal) = (d.ρ0,d.rs,d.α)
scale_radius(d::GeneralIsothermal) = d.rs

#### Evaluation

function ρ(d::GeneralIsothermal,r::Real)
    ρ0,rs,α = params(d)
    ρ0 * (r/rs)^-α
end
function invρ(d::GeneralIsothermal,x::Real)
    ρ0,rs,α = params(d)
    rs * (x/ρ0)^(-1/α)
end
function dρ_dr(d::GeneralIsothermal,r::Real)
    ρ0,rs,α = params(d)
    -α * ρ0 * (r/rs)^-(α+1) / rs
end
function Σ(d::GeneralIsothermal,r::Real)
    # This is the abel integral. Because it integrates to infinite along the line of sight, the integral of Σ,
    # e.g. quadgk(x->2π*x*Σ(d,x),0,R) will not equal 10 at Rmax, if you used the total mass enclosed within
    # Rmax to define the GeneralIsothermal via the GeneralIsothermal(rs::Real,α::Real,M::Real,Rmax::Real) method. 
    ρ0,rs,α = params(d)
    (α >= 2 || α <= 1) ? throw(DomainError(α,"Σ is only finite for 1<α<2.")) : sqrt(π) * r^(-α+1) * ρ0 * rs^α *
        gamma((α-1)/2) / gamma(α/2)
end
function invΣ(d::GeneralIsothermal,x::Real)
    ρ0,rs,α = params(d)
    if (α >= 2 || α <= 1)
        throw(DomainError(α,"Σ is only finite for 1<α<2."))
    elseif x<=0
        throw(DomainError(x,"invΣ only valid for x>0."))
    else
        (sqrt(π) * ρ0 * rs^α * gamma((α-1)/2) / gamma(α/2) / x)^(α-1)
    end
end
function dΣ_dr(d::GeneralIsothermal,r::Real)
    ρ0,rs,α = params(d)
    (α >= 2 || α <= 1) ? throw(DomainError(α,"Σ is only finite for 1<α<2.")) : -sqrt(π) * ρ0 * rs^α * (α-1) *
        gamma((α-1)/2) / gamma(α/2) / r^α
end
function M(d::GeneralIsothermal,r::Real)
    ρ0,rs,α=params(d)
    α >= 3 ? throw(DomainError(α,"Enclosed mass is only finite for α<3.")) : 4π * ρ0 * rs^α * r^(3-α) / (3 - α)
end
function dM_dr(d::GeneralIsothermal,r::Real)
    ρ0,rs,α=params(d)
    α >= 3 ? throw(DomainError(α,"Enclosed mass is only finite for α<3.")) : 4π * ρ0 * rs^α * r^(2-α)
end
# function Mtot(d::ExponentialDisk)
#     Σ0,rs=params(d)
#     2π * Σ0 * rs^2
# end
# function ccdf(d::ExponentialDisk,r::Real)
#     Σ0,rs=params(d)
#     ee = exp(-r/rs)
#     r*ee/rs + ee
# end
# cdf(d::ExponentialDisk,r::Real) = 1-ccdf(d,r)
# function cquantile(d::ExponentialDisk,x::Real)
#     Σ0,rs=params(d)
#     -rs * (1+lambertw(-x/ℯ,-1))
# end
# quantile(d::ExponentialDisk,x::Real) = cquantile(d,1-x)
