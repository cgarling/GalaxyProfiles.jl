# """
#     NFW(ρ0::Real, rs::Real)
#     NFW(ρ0::Unitful.Density, rs::Unitful.Length)

# Type describing the [Navarro-Frenk-White 1996 (NFW)](http://adsabs.harvard.edu/abs/1996ApJ...462..563N) density profile with scale radius `rs` and characteristic density `ρ0`. The density profile is

# ```math
# \\rho(r) = \\frac{\\rho_0}{(r/R_s) \\, (1+r/R_s)^2}
# ```

# The fields of `NFW` are `ρ0, rs`. The default units of `NFW` are `[ρ0] = [Msun/kpc^3], [r, rs] = [kpc], [M] = [Msun]`. This is important for quantities like [`Vcirc`](@ref), [`Vesc`](@ref), [`Φ`](@ref), [`∇Φ`](@ref), and [`∇∇Φ`](@ref) which involve the gravitational constant `G`; these will give incorrect results if the fields of `NFW` or the provided `r` are in different units.

# The following public methods are defined on this type:
#  - [`ρ`](@ref), [`invρ`](@ref), [`∇ρ`](@ref), [`ρmean`](@ref), [`Σ`](@ref), [`M`](@ref), [`∇M`](@ref), [`invM`](@ref), [`Mproj`](@ref), [`∇Mproj`](@ref), [`Vmax`](@ref), [`Vesc`](@ref), [`Φ`](@ref), [`∇Φ`](@ref), [`∇∇Φ`](@ref)
# """
# struct NFW{T<:Real} <: AbstractDensity
#     ρ0::T
#     rs::T
# end
# NFW(ρ0::Real, rs::Real) = NFW(promote(ρ0,rs)...)

"""
    CoreNFW(ρ0::Real, rs::Real, rc::Real, n::Real)

Type implementing the CoreNFW density profile of [Read et al. 2016](https://ui.adsabs.harvard.edu/abs/2016MNRAS.459.2573R/abstract). This is a modified [`NFW`](@ref) profile with lower central densities calibrated to reproduce the density profiles observed in their simulations. The parameters `ρ0` and `rs` are identical to those for the NFW profile. However, the CoreNFW profile makes the following alteration to the enclosed mass,

```math
\\begin{aligned}
M_\\text{c,NFW}(<R) &= M_\\text{NFW}(<R) \\times f^n \\newline
f^n &= \\left[ \\text{tanh} \\left( \\frac{r}{r_c} \\right) \\right]^n
\\end{aligned}
```

where the core radius ``r_c`` and ``n``, which controls how shallow the core is (no core with `n==0` and complete core with `n==1`), are new free parameters.

# See also
 - Constructor that uses galaxy properties, [`CoreNFWGalaxy`](@ref).
"""
struct CoreNFW{TNFW, Trc, Tn} <: AbstractDensity
    # We'll save the NFW object rather than ρ0 and rs so we don't have to create it
    # every time we want to evaluate the non-cored NFW profile
    NFW::TNFW
    rc::Trc
    n::Tn
end
function CoreNFW(ρ0, rs, rc, n)
    return CoreNFW(NFW(ρ0, rs), rc, n)
end
"""
    CoreNFWGalaxy(ρ0::Real, rs::Real, t_sf::Real, rhalf::Real;
                  κ::Real=4//100, η::Real=175//100)

An alternative parameterization for the [`CoreNFW`](@ref) profile based on galaxy properties.

The standard [`CoreNFW`](@ref) type uses ``n`` and ``r_c`` to parameterize the Read et al. 2016 density model, but the authors relate these halo-level parameters to galaxy-level parameters in the following ways. They find ``n`` to be related to the total star formation time,

```math
\\begin{aligned}
n &= \\text{tanh} \\left( q \\right) \\newline
q &= \\kappa \\, \\frac{t_\\text{SF}}{t_\\text{dyn}}
\\end{aligned}
```

where ``t_\\text{dyn}`` is the circular orbit time at the [`NFW`](@ref) profile scale radius ``r_s``. This is the orbital circumference divided by the speed; for profile `p`, this is `2π * p.rs / GalaxyProfiles.Vcirc(p, p.rs)`.

They also find the core radius ``r_c`` to be related to the projected stellar half-mass radius as

```math
r_c = \\eta \\, R_{h}
```

Thus if the total star formation time and projected stellar half-mass radius are known for a galaxy, one can calculate ``n`` and ``r_c`` given the alternative parameters ``\\eta`` and ``\\kappa``. In their simulations, Read et al. find ``\\kappa=0.04`` and ``\\eta=1.75``, so these are the default values for this alternate parameterization.
"""
function CoreNFWGalaxy(ρ0, rs, t_sf, rhalf;
                       κ=4//100, η=175//100)
    p = NFW(ρ0, rs)
    # Circular orbital time at NFW scale radius; is this just 4x our `dynamical_time`?
    t_dyn = rs * 2 * π / GalaxyProfiles.Vcirc(p, rs)
    # Convert from kpc / (km / s) to yr
    t_dyn *= 8202315653129241//8388608
    q = κ * t_sf / t_dyn # Equation 18
    n = tanh(q)          # Equation 18
    rc = η * rhalf       # Equation 20
    return CoreNFW(p, rc, n)
end


#### Parameters

params(d::CoreNFW) = (d.NFW.ρ0, d.NFW.rs, d.rc, d.n)
scale_radius(d::CoreNFW) = d.NFW.rs

#### Convenience functions

# Doesn't actually shorten anything 
# f_cnfw(r,rc) = tanh(r/rc) # Equation 17, Read 2016

#### Evaluation

function ρ(d::CoreNFW, r)
    rc, n = d.rc, d.n
    f = tanh(r/rc) # Equation 17, Read 2016
    fn = f^n
    prof = d.NFW
    # Equation 21,
    return ρ(prof, r) * fn + (n * (fn/f) * (1 - f^2)) / (rc * r^2 * 4 * π) * M(prof, r)
end
# function ρ(d::CoreNFW, r::Real)
#     ρ0, rs, rc, n = params(d)
#     x = r / rs
#     ρnfw = ρ0 / x / (1 + x)^2 # Standard NFW density
#     # Mnfw = ρ0 * 4 * π * rs^3 * NFWmu(r, rs)
#     f = tanh(r / rc) # Equation 17, Read 2016
#     fn = f^n
#     # prof = d.NFW
#     # Equation 21,
#     # return ρ(prof, r) * fn + (n * (fn/f) * (1 - f^2)) / (rc * r^2 * 4 * π) * M(prof, r)
#     return ρnfw * fn + (n * (fn/f) * (1 - f^2)) / (rc * r^2 * 4 * π) * M(prof, r)
# end
# function invρ(d::NFW, x::Real)
#     ρ0, rs = params(d)
#     tmp = cbrt(2x / (2x + 27ρ0 + 3*sqrt(3ρ0*(27ρ0+4x)))) # cbrt more efficient than ^(1/3)
#     rs/3 * (tmp + inv(tmp) - 2) 
# end
function ∇ρ(d::CoreNFW, r::Real)
    ρ0, rs, rc, n = params(d)
    return rs^3 * ρ0 / (r^3 * rc^2 * (r + rs)^3) * tanh(r/rc)^n *
        (-r * rc^2 * (3r + rs) + 8n * (r + rs) / (cosh(4r/rc) - 1) * 
        (r * (r + rs) * (n - cosh(2r/rc)) * (-r + (r + rs) * log(r/rs + 1)) -
        rc * (-r * (2r + rs) + (r + rs)^2 * log(r/rs + 1)) * sinh(2r/rc)))
end
# function ρmean(d::NFW, r::Real)
#     ρ0, rs = params(d)
#     3 * rs^3 * ρ0 * (-r + ( r + rs ) * log( ( r + rs ) /
#             rs ) ) / (r^3 * ( r + rs ) )
# end
# # invρmean fallback to common.jl is fine
# function Σ(d::NFW{T}, r::S) where {T, S<:Real}
#     U = promote_type(T, S)
#     ρ0, rs = params(d)
#     x = r/rs
#     if x ≈ 1 # approx for very near 1
#         U(2 * ρ0 * rs / 3)
#     elseif x < 1
#         x2m1 = x^2 - 1
#         rs * ρ0 * 2 / x2m1 * (1 - 2 / sqrt(-x2m1) * atanh( sqrt( (1 - x) / (x + 1) ) ) )
#     else
#         x2m1 = x^2 - 1
#         rs * ρ0 * 2 / x2m1 * (1 - 2 / sqrt(x2m1) * atan( sqrt( (x - 1) / (x + 1) ) ) )
#     end
# end
# # ∇Σ
# # Σmean fallback to common.jl is fine
# # invΣ
function M(d::CoreNFW, r::Real)
    rc, n = d.rc, d.n
    f = tanh(r/rc) # Equation 17, Read 2016
    fn = f^n
    return M(d.NFW, r) * fn # Equation 16
end
function ∇M(d::CoreNFW, r::Real)
    rc, n = d.rc, d.n
    x = r/rc
    f = tanh(x) # Equation 17, Read 2016
    fn = f^n
    prof = d.NFW
    # Product rule
    return fn * ∇M(prof, r) + n * sech(x)^2 * (fn/f) / rc * M(prof, r)
end
# function invM(d::NFW{T}, x::S) where {T, S<:Real}
#     U = promote_type(T, S)
#     ρ0,rs = params(d)
#     ρ0, rs, x = promote(ρ0, rs, x)
#     t = -1/lambertw(-exp(-1-(x/(4 * U(π) * rs^3 * ρ0))))
#     (t - 1) * rs
# end
# # Mtot
# function Mproj(d::NFW{T}, r::S) where {T, S<:Real}
#     U = promote_type(T, S)
#     ρ0, rs = params(d)
#     ρ0, rs, r = promote(ρ0, rs, r)
#     x = r/rs
#     if x ≈ 1
#         4 * U(π) * ρ0 * rs^3 * ( x + log( x / 2 ) )
#     elseif x<1
#         return 4 * U(π) * ρ0 * rs^3 * ( acosh( 1 / x ) / sqrt( 1 - x^2 ) + log( x / 2 ) )
#     else
#         return 4 * U(π) * ρ0 * rs^3 * ( acos( 1 / x ) / sqrt( x^2 - 1 ) + log( x / 2 ) )
#     end
# end
# function ∇Mproj(d::NFW{T},r::S) where {T,S<:Real}
#     U = promote_type(T,S)
#     ρ0,rs = params(d)
#     ρ0,rs,r = promote(ρ0,rs,r)
#     x = r/rs
#     if x ≈ 1 
#         4 * U(π) * ρ0 * rs^2 * (1+inv(x)) / 6
#     elseif x<1
#         ix = inv(x)
#         # return 4 * U(π) * ρ0 * rs^2 * (ix + sqrt(ix+ix^2)/(x-1)/(1+x)^(3/2) + x*asech(x)/(1-x^2)^(3/2))
#         return 4 * U(π) * ρ0 * rs^2 * (ix + sqrt(ix+ix^2)/(x-1)/sqrt((1+x)^3) + x*asech(x)/sqrt((1-x^2)^3)) # sqrt(x^3) is much more efficient than x^(3/2), annoyingly
#     else
#         ix = inv(x)
#         # return 4 * U(π) * ρ0 * rs^2 * ( ix + ix^2/sqrt(1-ix^2)/sqrt(x^2-1) - x*acos(ix)/(x^2-1)^(3/2))
#         return 4 * U(π) * ρ0 * rs^2 * ( ix + ix^2/sqrt(1-ix^2)/sqrt(x^2-1) - x*acos(ix)/sqrt((x^2-1)^3))
#     end
# end
# # invMproj fallback to common.jl is fine
# # Vcirc fallback to common.jl is fine
# # Vesc fallback to common.jl is fine
# function Vmax(d::NFW{T}) where T
#     ρ0,rs = params(d)
#     r = T(constants.nfwvmaxpar) * rs
#     return Vcirc(d,r), r
# end
# SatGen has an approximate potential at https://github.com/shergreen/SatGen/blob/7f76cf86ead6fc8e1e608e0353c2c8e05d4bff21/profiles.py#L1050
# function Φ(d::NFW{T}, r::S) where {T, S<:Real}
#     U = promote_type(T, S)
#     ρ0, rs = params(d)
#     ρ0, rs, r = promote(ρ0, rs, r)
#     4 * U(π) * U(constants.Gvelkpc) * ρ0 * rs^3 * (log(rs) - log(r+rs)) / r
# end
# function ∇Φ(d::NFW{T}, r::S) where {T, S<:Real}
#     U = promote_type(T, S)
#     ρ0, rs = params(d)
#     ρ0, rs, r = promote(ρ0, rs, r)
#     4 * U(π) * U(constants.Gvelkpc2) * rs^3 * ρ0 * (log(r+rs) - log(rs) - r/(r+rs)) / r^2
# end
# function ∇∇Φ(d::NFW{T}, r::S) where {T, S<:Real}
#     U = promote_type(T, S)
#     ρ0,rs = params(d)
#     ρ0, rs, r = promote(ρ0, rs, r)
#     4 * U(π) * U(constants.Gvelkpc2) * rs^3 * ρ0 * (2*log(rs) - 2*log(r+rs) + (3*r^2 + 2*rs*r)/(r+rs)^2) / r^3
# end
# # Could add kinetic energy, line-of-sight velocity dispersion, maybe a few other quantities that could be useful from python/profiles/nfw.py
