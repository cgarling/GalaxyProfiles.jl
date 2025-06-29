"""
    Plummer(M::Real, a::Real)
    Plummer(M::Unitful.Mass, a::Unitful.Length)

Type implementing the density profile first proposed by [Plummer 1911](http://adsabs.harvard.edu/abs/1911MNRAS..71..460P) (see also [Dejonghe 1987](https://ui.adsabs.harvard.edu/abs/1987MNRAS.224...13D/abstract)); this density profile is defined by the potential

```math
\\Phi(r) = -\\frac{G \\, M}{\\sqrt{ r^2 + a^2 } }
```

which corresponds to a density profile of

```math
\\rho(r) = \\frac{3 \\, M}{4 \\pi a^3} \\, \\left( 1 + \\frac{r^2}{a^2} \\right)^{-5/2}
```

where `M` is the total mass of the system and `a` is the characteristic scale radius of the system. Following this, the fields of `Plummer` are `M, a`. This scale radius can be converted to the 3D half-light (or half-mass) radius via

```math
r_h = \\frac{a}{\\sqrt{ 0.5^{-2/3} - 1} }
```

the method [`GalaxyProfiles.plummer_a_to_rh`](@ref) is provided to perform this conversion. In projection (i.e., along a line-of-sight) the half-light radius is equal to the Plummer scale radius `a`, verifiable by [`quantile2D(d::Plummer, 0.5) == scale_radius(d)`](@ref quantile2D). 

The default units of `Plummer` are `[M] = [Msun], [a, r] = [kpc]`. This is important for quantities like [`Vcirc`](@ref), [`Vesc`](@ref), [`Φ`](@ref), [`∇Φ`](@ref), and [`∇∇Φ`](@ref) which involve the gravitational constant `G`; these will give incorrect results if the fields of `Plummer` or the provided `r` are in different units.

The following public methods are defined on this type:
 - [`Mtot`](@ref), [`ρ`](@ref), [`invρ`](@ref), [`∇ρ`](@ref), [`ρmean`](@ref), [`invρmean`](@ref), [`Σ`](@ref), [`∇Σ`](@ref), [`Σmean`](@ref), [`invΣ`](@ref), [`M`](@ref), [`∇M`](@ref), [`invM`](@ref), [`Mproj`](@ref), [`∇Mproj`](@ref), [`invMproj`](@ref), [`Vcirc`](@ref), [`Vesc`](@ref), [`Φ`](@ref), [`∇Φ`](@ref), [`∇∇Φ`](@ref), [`cdf2D`](@ref), [`cdf3D`](@ref), [`ccdf2D`](@ref), [`ccdf3D`](@ref), [`quantile2D`](@ref), [`quantile3D`](@ref), [`cquantile2D`](@ref), [`cquantile3D`](@ref).
"""
struct Plummer{T <: Real} <: AbstractDensity{T}
    M::T
    a::T
end
Plummer(M::Real, a::Real) = Plummer(promote(M,a)...)

#### Parameters

params(d::Plummer) = (d.M, d.a)
scale_radius(d::Plummer) = d.a
Mtot(d::Plummer) = d.M

#### Convenience functions
"""
    plummer_a_to_rh(a::Real)

Returns the 3-D half-light (or half-mass) radius given the Plummer scale radius `a`. This is equivalent to [`quantile3D(d::Plummer, 0.5)`](@ref quantile3D) but faster because the argument `x=0.5` is known at compile-time. Note that in projection (i.e., along a line-of-sight) the half-light radius is simply equal to the Plummer scale radius `a`, verifiable by [`quantile2D(d::Plummer, 0.5) == scale_radius(d)`](@ref quantile2D).
"""
plummer_a_to_rh(a::T) where T <: Real = a / sqrt(cbrt(convert(float(T), 1//2))^-2 - 1) # The constant denominator is ~inv(1.3).
"""
    plummer_rh_to_a(a::Real)

Returns the Plummer scale radius `a` given a 3-D half-light (or half-mass) radius. Note that in projection (i.e., along a line-of-sight) the half-light radius is simply equal to the Plummer scale radius `a`, verifiable by [`quantile2D(d::Plummer, 0.5) == scale_radius(d)`](@ref quantile2D).
"""
plummer_rh_to_a(rh::T) where T <: Real = rh * sqrt(cbrt(convert(float(T), 1//2))^-2 - 1)
"""
    plummer_angular_avalue(absmag, sb, DM)

This is an observational utility designed for use with the astronomical magnitude system. Returns the Plummer scale radius `a` given a desired absolute magnitude `absmag`, average surface brightness within the half-light radius `sb`, and distance modulus `DM`. If the units of `sb` are, for example, [mags/arcsec^2], and the units of `absmag` and `DM` should always be magnitudes, then this will return the Plummer scale radius `a` in [arcsec]. Mathematically, we solve the equation for the surface brightness `S` as a function of magnitude `m = absmag + DM` and area `A=πr^2`, with `r` being the half-light radius:

```math
\\begin{aligned}
S &= m + 2.5 \\times \\log(A) = m + 2.5 \\times \\log(πr^2) = m + 2.5 \\times ( \\log(π) + 2\\log(r) ) \\newline
r &= \\text{exp10} \\left(  \\frac{S - m - 2.5 \\, \\log(π)}{5} \\right) = \\text{exp10} \\left(  \\frac{S - m}{5} \\right) \\times π^{1/2}
\\end{aligned}
```

We then only need to convert the half-light radius `r` to the scale radius `a`, which is just the inverse of [`GalaxyProfiles.plummer_a_to_rh`](@ref).

# Examples
For a Plummer profile with a total magnitude of -5 (`flux = -2.5*log10(mag) = 100`), an average surface brightness within the half-light radius of 25 mag/arcsec^2, and a distance of 1 Mpc (`DM = 5*log10(1e6 [pc] / 10 [pc]) = 25`), we can compute the scale radius of the corresponding Plummer profile in arcseconds as 
```jldoctest
julia> isapprox( GalaxyProfiles.plummer_angular_avalue(-5, 25, 25), 4.32; rtol=1e-2)
true
```
"""
function plummer_angular_avalue(absmag::T, sb::S, DM::V) where {T,S,V}
    U = float(promote_type(T,S,V)) # Must be a float type
    return exp10((sb - (absmag + DM)) / 5) / sqrt(U(π)) * sqrt(cbrt(U(1//2))^-2  - 1)
end
"""
    plummer_unscaled_density(r, M, a)

Returns the unscaled density for a Plummer profile at radius `r` with total mass `M` and Plummer scale radius `a`. 
"""
plummer_unscaled_density(r, M, a) = @fastpow sqrt((1+(r/a)^2)^-5)
"""
    plummer_unscaled_density_deriv(r, M, a)

Returns the unscaled radial derivative of the density for a Plummer profile at radius `r` with total mass `M` and Plummer scale radius `a`. 
"""
plummer_unscaled_density_deriv(r, M, a) = @fastpow -5r / (a^2 * sqrt((1+(r/a)^2)^7))
plummer_unscaled_ρmean(r, M, a) = @fastpow a^4 * sqrt(1 + (r/a)^2) / (a^2 + r^2)^2
plummer_unscaled_Σ(r, M, a) = @fastpow 4a^5 / 3(a^2 + r^2)^2
plummer_unscaled_∇Σ(r, M, a) = @fastpow -16a^5 * r / 3(a^2 + r^2)^3
# plummer_unscaled_Σmean(r, M, a) = 2a^2 * (1 - (1 + (r/a)^2)^(-3/2)) / (3r^2)
plummer_unscaled_Σmean(r, M, a) = 4a^3 / 3(a^2 + r^2)
# Dont really think these are that useful. Maybe just remove later.

#### Evaluation

# function ρ(d::Plummer{T}, r::S) where {T,S} # = (3*self.M/(4.*np.pi*self.a**3.)) * self.unscaled_density(r)
function ρ(d::Plummer, r::Real) # = (3*self.M/(4.*np.pi*self.a**3.)) * self.unscaled_density(r)
    r, M, a = promote(r, Mtot(d), scale_radius(d))
    return 3M / (a^3 * 4 * π) * plummer_unscaled_density(r, M, a)
end
function invρ(d::Plummer{T}, x::S) where {T <: Real, S <: Real}
    x, M, a = promote(x, Mtot(d), scale_radius(d))
    U = float(promote_type(T,S))
    return a * sqrt( U(2^(-4/5) * 3^(2/5) * π^(-2/5)) * (a^3 * x / M)^(-2//5) - 1 )
end
function ∇ρ(d::Plummer, r::Real)
    r, M, a = promote(r, Mtot(d), scale_radius(d))
    return 3M / (a^3 * 4 * π) * plummer_unscaled_density_deriv(r, M, a)
end
function ρmean(d::Plummer, r::Real)
    r, M, a = promote(r, Mtot(d), scale_radius(d))
    return 3M / (a^3 * 4 * π) * plummer_unscaled_ρmean(r, M, a)
end
function invρmean(d::Plummer, x::Real)
    x, M, a = promote(x, Mtot(d), scale_radius(d))
    return sqrt((M * 6 / (π * x))^(2//3) - 4a^2) / 2
end
function Σ(d::Plummer, r::Real)
    r, M, a = promote(r, Mtot(d), scale_radius(d))
    return M * a^2 / π / (a^2 + r^2)^2
end
function ∇Σ(d::Plummer, r::Real)
    r, M, a = promote(r, Mtot(d), scale_radius(d))
    return -a^2 * M * r * 4 / π / (a^2 + r^2)^3
end
function Σmean(d::Plummer, r::Real)
    r, M, a = promote(r, Mtot(d), scale_radius(d))
    return M / (a^2 + r^2) / π
end
function invΣ(d::Plummer, x::Real)
    x, M, a = promote(x, Mtot(d), scale_radius(d))
    return a * sqrt(sqrt(M / x / π) / a - 1)
end
# function M(d::Plummer, r::Real)
#     r, M, a = promote(r, Mtot(d), scale_radius(d))
#     return a * M * r^3 * sqrt(1 + (r/a)^2) / (a^2 + r^2)^2
# end
function ∇M(d::Plummer, r::Real)
    r, M, a = promote(r, Mtot(d), scale_radius(d))
    return 3 * a * M * r^2 / (a^2 + r^2)^2 / sqrt(1 + (r/a)^2)
end
function invM(d::Plummer, x::Real) # Actually basing this off quantile3D from Aarseth 1974.
    x, M, a = promote(x, Mtot(d), scale_radius(d))
    return scale_radius(d) / sqrt(cbrt((x/M)^-2) - 1)
end
function Mproj(d::Plummer, r::Real)
    r, M, a = promote(r, Mtot(d), scale_radius(d))
    return M * r^2 / (a^2 + r^2)
end
function ∇Mproj(d::Plummer, r::Real)
    M, a = Mtot(d), scale_radius(d)
    return 2a^2 * M * r / (a^2 + r^2)^2
end
function invMproj(d::Plummer, x::Real)
    M, a = Mtot(d), scale_radius(d)
    return a * sqrt(x) / sqrt(M - x)
end
# Vcirc fallback to common.jl is fine
# Vesc fallback to common.jl is fine
# There is a variation of this for all hypervirial density-potential pairs; see Evans2005
function σr(d::Plummer, r::T, β::S) where {T <: Real, S <: Real}
    U = promote_type(T, S)
    # HypergeometricFunctions._₂F₁ will return NaN if z (last argument) is Inf (r=0), but it should return 1...
    int_factor = r ≈ zero(T) ? one(U) : _₂F₁(-β, one(U), 4-β, -r^-2)
    return sqrt(-Φ(d, r) / (6-2β) * int_factor)
end
function Φ(d::Plummer{T}, r::S) where {T, S <: Real}
    r, M, a = promote(r, Mtot(d), scale_radius(d))
    U = float(promote_type(T,S))
    return -U(constants.Gvelkpc) * M / sqrt(r^2 + a^2)
end
function ∇Φ(d::Plummer{T}, r::S) where {T, S <: Real}
    r, M, a = promote(r, Mtot(d), scale_radius(d))
    U = float(promote_type(T,S))
    return U(constants.Gvelkpc2) * M * r / sqrt((r^2 + a^2)^3)
end
function ∇∇Φ(d::Plummer{T}, r::S) where {T, S <: Real}
    r, M, a = promote(r, Mtot(d), scale_radius(d))
    U = float(promote_type(T,S))
    denom = a^2 + r^2
    denom2 = denom^2
    return U(constants.Gvelkpc2) * M * (a^2 - 2r^2) / sqrt(denom2 * denom2 * denom)
end
cdf2D(d::Plummer, r::Real) = r^2 / (scale_radius(d)^2 + r^2) # Just Mproj(d,R) / Mtot(d). 
# ccdf2D fallback to 1 - cdf2D(d,r) in common.jl is fine. 
function cdf3D(d::Plummer, r::Real) # Just M(d,R) / Mtot(d). 
    a = scale_radius(d)
    return a * r^3 * sqrt(1 + (r/a)^2) / (a^2 + r^2)^2
end
# ccdf3D fallback to 1 - cdf3D(d,r) in common.jl is fine. 
quantile2D(d::Plummer, x::Real) = scale_radius(d) * sqrt(x) / sqrt(1-x)
cquantile2D(d::Plummer, x::Real) = quantile2D(d, 1-x)
quantile3D(d::Plummer, x::Real) = scale_radius(d) / sqrt(cbrt(x^-2) - 1) # This is the Aarseth 1974 result
# function quantile3D(d::Plummer, x::Real) # This is the full expansion but it is virtually the same as above. 
#     a = scale_radius(d)
#     a2 = a^2
#     a4 = a2^2
#     a8 = a4^2
#     a12 = a8 * a4
#     x2 = x^2
#     x4 = x2^2
#     x6 = x4 * x2
#     return sqrt( ( cbrt(2) * a4 * x2 / ( (x2 - 1) * cbrt( sqrt(a12 * x4 - 2 * a12 * x6 + a12 * x4^2) - a4 * a2 * x2 - a4 * a2 * x4) ) ) +
#         ( cbrt( -a4 * a2 * x2 - a4 * a2 * x4 + sqrt(a12 * x4 - 2 * a12 * x6 + a12 * x4^2) ) / ( cbrt(2) * (x2 - 1) ) ) -
#         a2 * x2 / (x2 - 1) )
# end
cquantile3D(d::Plummer, x::Real) = quantile3D(d, 1-x)
