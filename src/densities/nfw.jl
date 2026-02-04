"""
    NFW(ρ0::Real, rs::Real)
    NFW(ρ0::Unitful.Density, rs::Unitful.Length)

Type describing the [Navarro-Frenk-White 1996 (NFW)](http://adsabs.harvard.edu/abs/1996ApJ...462..563N) density profile with scale radius `rs` and characteristic density `ρ0`. The density profile is

```math
\\rho(r) = \\frac{\\rho_0}{(r/R_s) \\, (1+r/R_s)^2}
```

The fields of `NFW` are `ρ0, rs`. The default units of `NFW` are `[ρ0] = [Msun/kpc^3], [r, rs] = [kpc], [M] = [Msun]`. This is important for quantities like [`Vcirc`](@ref), [`Vesc`](@ref), [`Φ`](@ref), [`∇Φ`](@ref), and [`∇∇Φ`](@ref) which involve the gravitational constant `G`; these will give incorrect results if the fields of `NFW` or the provided `r` are in different units.

The following public methods are defined on this type:
 - [`ρ`](@ref), [`invρ`](@ref), [`∇ρ`](@ref), [`ρmean`](@ref), [`Σ`](@ref), [`M`](@ref), [`∇M`](@ref), [`invM`](@ref), [`Mproj`](@ref), [`∇Mproj`](@ref), [`Vmax`](@ref), [`Vesc`](@ref), [`Φ`](@ref), [`∇Φ`](@ref), [`∇∇Φ`](@ref)
"""
struct NFW{Tρ, Tr} <: AbstractDensity
    ρ0::Tρ
    rs::Tr
end

#### Parameters

params(d::NFW) = (d.ρ0, d.rs)
scale_radius(d::NFW) = d.rs

#### Convenience functions

NFWmu(x) = log(1+x) - x / (1+x)
NFWmu(r, rs) = NFWmu(r/rs)

#### Evaluation

function ρ(d::NFW, r)
    ρ0,rs = params(d)
    x = r/rs
    return ρ0 / x / (1+x)^2
end
function invρ(d::NFW, x::Real)
    ρ0, rs = params(d)
    tmp = cbrt(2x / (2x + 27ρ0 + 3sqrt(3ρ0*(27ρ0+4x)))) # cbrt more efficient than ^(1/3)
    return rs/3 * (tmp + inv(tmp) - 2) 
end
function ∇ρ(d::NFW, r::Real)
    ρ0, rs = params(d)
    return -ρ0 * rs^3 * (3r * rs) / (r^2 * (r + rs)^3)
end
function ρmean(d::NFW, r::Real)
    ρ0, rs = params(d)
    return 3rs^3 * ρ0 * (-r + (r+rs) * log((r+rs) / rs)) /
        (r^3 * (r+rs))
end
# invρmean fallback to common.jl is fine
function Σ(d::NFW{T}, r::S) where {T, S}
    U = promote_type(T, S)
    ρ0, rs = params(d)
    x = r/rs
    if x ≈ 1 # approx for very near 1
        return U(2 * ρ0 * rs / 3)
    elseif x < 1
        x2m1 = x^2 - 1
        return rs * ρ0 * 2 / x2m1 * (1 - 2 / sqrt(-x2m1) * atanh( sqrt((1-x) / (x+1)) ) )
    else
        x2m1 = x^2 - 1
        return rs * ρ0 * 2 / x2m1 * (1 - 2 / sqrt(x2m1) * atan( sqrt((x - 1) / (x + 1)) ) )
    end
end
# ∇Σ
# Σmean fallback to common.jl is fine
# invΣ
function M(d::NFW, r::Real)
    ρ0, rs = params(d)
    return rs^3 * ρ0 * NFWmu(r, rs) * 4 * π
end
function ∇M(d::NFW, r::Real)
    ρ0, rs = params(d)
    return r * rs * ρ0 * 4 * π / (1 + r/rs)^2
end
function invM(d::NFW, x::Real)
    ρ0, rs = params(d)
    t = -inv(lambertw(-exp(-1-(x/(rs^3 * ρ0 * 4 * π)))))
    return (t - 1) * rs
end
# Mtot
function Mproj(d::NFW, r::Real)
    ρ0, rs = params(d)
    x = r/rs
    if x ≈ 1
        return ρ0 * rs^3 * 4 * π * (x + log(x/2))
    elseif x<1
        return ρ0 * rs^3 * 4 * π * (acosh(1/x) / sqrt(1-x^2) + log(x/2))
    else
        return ρ0 * rs^3 * 4 * π * (acos(1/x) / sqrt(x^2-1) + log(x/2))
    end
end
function ∇Mproj(d::NFW, r::Real)
    ρ0, rs = params(d)
    x = r/rs
    invx = inv(x)
    if x ≈ 1 
        return ρ0 * rs^2 * 4 * π * (1+invx) / 6
    elseif x < 1
        # return 4 * U(π) * ρ0 * rs^2 * (ix + sqrt(ix+ix^2)/(x-1)/(1+x)^(3/2) + x*asech(x)/(1-x^2)^(3/2))
        return ρ0 * rs^2 * 4 * π * (invx + sqrt(invx+invx^2)/(x-1)/sqrt((1+x)^3) + x*asech(x)/sqrt((1-x^2)^3)) # sqrt(x^3) is much more efficient than x^(3/2), annoyingly
    else
        # return 4 * U(π) * ρ0 * rs^2 * ( ix + ix^2/sqrt(1-ix^2)/sqrt(x^2-1) - x*acos(ix)/(x^2-1)^(3/2))
        return ρ0 * rs^2 * 4 * π * (invx + invx^2/sqrt(1-invx^2)/sqrt(x^2-1) - x*acos(invx)/sqrt((x^2-1)^3))
    end
end
# invMproj fallback to common.jl is fine
# Vcirc fallback to common.jl is fine
# Vesc fallback to common.jl is fine
function Vmax(d::NFW{T}) where T
    ρ0, rs = params(d)
    r = T(constants.nfwvmaxpar) * rs
    return Vcirc(d,r), r
end
# Roughly 2x faster than generic fallback
function σr(d::NFW{T}, r::S, β::S) where {T, S}
    U = promote_type(T, S)
    ρ0, rs = params(d)
    x = r / rs
    Φ0 = ρ0 * rs^2 * 4 * π * U(constants.Gvelkpc)
    return sqrt(Φ0 * quadgk(x->x^(2β-3) * NFWmu(x) / (1+x)^2, x, utilities.get_inf(x))[1] /
        x^(2β-1) * (1+x)^2)
end
function Φ(d::NFW{T}, r::S) where {T, S}
    U = promote_type(T, S)
    ρ0, rs = params(d)
    return U(constants.Gvelkpc) * ρ0 * rs^3 * 4 * π * (log(rs) - log(r+rs)) / r
end
function ∇Φ(d::NFW{T}, r::S) where {T, S}
    U = promote_type(T, S)
    ρ0, rs = params(d)
    return U(constants.Gvelkpc2) * rs^3 * ρ0 * 4 * π * (log(r+rs) - log(rs) - r/(r+rs)) / r^2
end
function ∇∇Φ(d::NFW{T}, r::S) where {T, S}
    U = promote_type(T, S)
    ρ0, rs = params(d)
    return U(constants.Gvelkpc2) * rs^3 * ρ0 * 4 * π * (2*log(rs) - 2*log(r+rs) + (3*r^2 + 2*rs*r)/(r+rs)^2) / r^3
end
# Could add kinetic energy, line-of-sight velocity dispersion, maybe a few other quantities that could be useful from python/profiles/nfw.py
