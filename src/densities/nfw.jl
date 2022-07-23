"""
    NFW(ρ0::Real,rs::Real)

Type describing the Navarro-Frenk-White (NFW) density profile with scale radius `rs` and characteristic density `ρ0`. The surface density profile is

```math
\\rho(r) = \\frac{\\rho_0}{(r/R_s) \\, (1+r/R_s)^2}
```

The fields of `NFW` are `ρ0, rs`. The default units of `NFW` are `[ρ0] = [Msun/kpc^3], [r, rs] = [kpc], [M] = [Msun]`. This is important for quantities like `Vcirc`, `Φ`, `∇Φ` which involve `G`; these will give incorrect results if the fields of `NFW` or the provided `r` are in different units.
"""
struct NFW{T<:Real} <: AbstractDensity
    ρ0::T
    rs::T
end
NFW(ρ0::Real,rs::Real) = NFW(promote(ρ0,rs)...)

#### Parameters

params(d::NFW) = (d.ρ0,d.rs)
scale_radius(d::NFW) = d.rs

#### Convenience functions

NFWmu(r,rs) = (x=r/rs; log(1+x) - x / (1+x))
NFWmu(x) = log(1+x) - x / (1+x)

#### Evaluation

function ρ(d::NFW,r::Real)
    ρ0,rs = params(d)
    x = r/rs
    ρ0 / x / (1+x)^2
end
function invρ(d::NFW,x::Real)
    ρ0,rs = params(d)
    tmp = (2x / (2x + 27ρ0 + 3*sqrt(3ρ0*(27ρ0+4x))))^(1/3)
    rs/3 * (tmp + inv(tmp) - 2) 
end
function ∇ρ(d::NFW,r::Real)
    ρ0,rs = params(d)
    -ρ0 * rs^3 * (3 * r * rs) / (r^2 * (r + rs)^3)
end
function ρmean(d::NFW,r::Real)
    ρ0,rs = params(d)
    3 * rs^3 * ρ0 * (-r + ( r + rs ) * log( ( r + rs ) /
            rs ) ) / (r^3 * ( r + rs ) )
end
function M(d::NFW,r::Real)
    ρ0,rs = params(d)
    4π * rs^3 * ρ0 * NFWmu(r,rs)
end
function invM(d::NFW,x::Real)
    ρ0,rs = params(d)
    t = -1/lambertw(-exp(-1-(x/(4π * rs^3 * ρ0))))
    (t-1)*rs
end
function ∇M(d::NFW,r::Real)
    ρ0,rs = params(d)
    4π * r * rs * ρ0 / (1+r/rs)^2
end
# Mtot
function Σ(d::NFW,r::Real)
    ρ0,rs = params(d)
    x = r/rs
    if x ≈ 1 # approx for very near 1
        2 * ρ0 * rs / 3
    elseif x < 1
        x2m1 = x^2 - 1
        rs * ρ0 * 2 / x2m1 * (1 - 2 / sqrt(-x2m1) * atanh( sqrt( (1 - x) / (x + 1) ) ) )
    else
        x2m1 = x^2 - 1
        rs * ρ0 * 2 / x2m1 * (1 - 2 / sqrt(x2m1) * atan( sqrt( (x - 1) / (x + 1) ) ) )
    end
end
# ∇Σ
# Σmean
# invΣ
function Mproj(d::NFW,r::Real)
    ρ0,rs = params(d)
    x = r/rs
    if x ≈ 1
        4π * ρ0 * rs^3 * ( x + log( x / 2 ) )
    elseif x<1
        return 4π * ρ0 * rs^3 * ( acosh( 1 / x ) / sqrt( 1 - x^2 ) + log( x / 2 ) )
    else
        return 4π * ρ0 * rs^3 * ( acos( 1 / x ) / sqrt( x^2 - 1 ) + log( x / 2 ) )
    end
end
function ∇Mproj(d::NFW,r::Real)
    ρ0,rs = params(d)
    x = r/rs
    if x ≈ 1 
        4π * ρ0 * rs^2 * (1+inv(x)) / 6
    elseif x<1
        ix = inv(x)
        return 4π * ρ0 * rs^2 * (ix + sqrt(ix+ix^2)/(x-1)/(1+x)^(3/2) + x*asech(x)/(1-x^2)^(3/2))
    else
        ix = inv(x)
        return 4π * ρ0 * rs^2 * ( ix + ix^2/sqrt(1-ix^2)/sqrt(x^2-1) - x*acos(ix)/(x^2-1)^(3/2))
    end
end
# invMproj
# need to think about the potential implementations and whether we want to use units or not
