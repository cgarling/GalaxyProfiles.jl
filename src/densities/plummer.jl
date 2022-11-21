"""
    Plummer(M::Real, a::Real)
    Plummer(M::Unitful.Mass, a::Unitful.Length)

Type implementing the density profile first proposed by [Plummer 1911](http://adsabs.harvard.edu/abs/1911MNRAS..71..460P); this density profile is defined by the potential

```math
\\Phi(r) = -\\frac{G \\, M}{\\sqrt{ r^2 + a^2 } }
```

which corresponds to a density profile of

```math
\\rho(r) = \\frac{3 \\, M}{4 \\pi a^3} \\, \\left( 1 + \\frac{r^2}{a^2} \\right)^{-5/2}
```

where `M` is the total mass of the system. Following this, the fields of `Plummer` are `M, a`, where `a` is the Plummer scale radius. This scale radius can be converted to the half-light (or half-mass) radius via

```math
r_h = \\frac{a}{\\sqrt{ 0.5^{-2/3} - 1} }
```

the method [`GalaxyProfiles.plummer_a_to_rh`](@ref) is provided to perform this conversion.

The default units of `Plummer` are `[M] = [Msun], [a, r] = [kpc]`. This is important for quantities like [`Vcirc`](@ref), [`Φ`](@ref), [`∇Φ`](@ref), and [`∇∇Φ`](@ref) which involve the gravitational constant `G`; these will give incorrect results if the fields of `Plummer` or the provided `r` are in different units.

The following methods are defined on this type:
 - [`Mtot`](@ref), [`ρ`](@ref), [`invρ`](@ref), [`∇ρ`](@ref), [`ρmean`](@ref), [`invρmean`](@ref), [`Σ`](@ref), [`∇Σ`](@ref), [`Σmean`](@ref), [`invΣ`](@ref), [`M`](@ref), [`∇M`](@ref), [`invM`](@ref), [`Mproj`](@ref), [`∇Mproj`](@ref), [`invMproj`](@ref), [`Vcirc`](@ref), [`Vesc`](@ref), [`Φ`](@ref), [`∇Φ`](@ref), [`∇∇Φ`](@ref)
"""
struct Plummer{T<:Real} <: AbstractDensity
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
    plummer_a_to_rh(a)

Returns the half-light (or half-mass) radius given the Plummer scale radius `a`.
"""
plummer_a_to_rh(a) = a / sqrt( inv( cbrt(0.5)^2 ) - 1 )
"""
    plummer_unscaled_density(r, M, a)

Returns the unscaled density for a Plummer profile at radius `r` with total mass `M` and Plummer scale radius `a`. 
"""
plummer_unscaled_density(r, M, a) =
    (interior = (1 + (r/a)^2); interior2 = interior^2; sqrt( inv( interior2^2 * interior ) ) )
    # (1+(r/a)^2)^(-5/2)
"""
    plummer_unscaled_density_deriv(r, M, a)

Returns the unscaled radial derivative of the density for a Plummer profile at radius `r` with total mass `M` and Plummer scale radius `a`. 
"""
plummer_unscaled_density_deriv(r, M, a) =
    (interior = (1 + (r/a)^2); interior2 = interior^2; interior4 = interior2^2;
     -5r / a^2 / sqrt( interior4 * interior2 * interior) )
# plummer_unscaled_density_deriv(r, M, a) = -5r / (a^2 * (1+(r/a)^2)^(7/2))
plummer_unscaled_ρmean(r, M, a) = a^4 * sqrt( 1 + (r/a)^2 ) / (a^2 + r^2)^2
plummer_unscaled_Σ(r, M, a) = (a2 = a^2; 4/3 * a2*a2*a / (a2 + r^2)^2)
# plummer_unscaled_Σ(r, M, a) = 4/3 * a^5 / (a^2 + r^2)^2
plummer_unscaled_∇Σ(r, M, a) = (a2 = a^2; -16/3 * a2 * a2 * a * r / (a2 + r^2)^3)
# plummer_unscaled_∇Σ(r, M, a) = -16/3 * a^5 * r / (a^2 + r^2)^3
# plummer_unscaled_Σmean(r, M, a) = 2a^2 * (1 - (1 + (r/a)^2)^(-3/2)) / (3r^2)
plummer_unscaled_Σmean(r, M, a) = 4/3 * a^3 / (a^2 + r^2)
# Dont really think these are that useful. Maybe just remove later.

#### Evaluation

function ρ(d::Plummer, r::Real) # = (3*self.M/(4.*np.pi*self.a**3.)) * self.unscaled_density(r)
    M, a = Mtot(d), scale_radius(d)
    return 3M / (4π*a^3) * plummer_unscaled_density(r, M, a)
end
function invρ(d::Plummer, x::Real)
    M, a = Mtot(d), scale_radius(d)
    return a * sqrt( 2^(-4/5) * 3^(2/5) * π^(-2/5) * (a^3 * x / M)^(-2/5) - 1 )
end
function ∇ρ(d::Plummer, r::Real)
    M, a = Mtot(d), scale_radius(d)
    return 3M / (4π*a^3) * plummer_unscaled_density_deriv(r, M, a)
end
function ρmean(d::Plummer, r::Real)
    M, a = Mtot(d), scale_radius(d)
    return 3M / (4π*a^3) * plummer_unscaled_ρmean(r, M, a)
end
function invρmean(d::Plummer, x::Real)
    M, a = Mtot(d), scale_radius(d)
    return sqrt( M^(2/3) * (6/π)^(2/3) / x^(2/3) - 4a^2 )
end
function Σ(d::Plummer, r::Real)
    M, a = Mtot(d), scale_radius(d)
    return M/π * a^2 / (a^2 + r^2)^2 # return 3M / (4π*a^3) * plummer_unscaled_Σ(r, M, a)
end
function ∇Σ(d::Plummer, r::Real)
    M, a = Mtot(d), scale_radius(d)
    return -4/π * a^2 * M * r / (a^2 + r^2)^3 # return 3M / (4π*a^3) * plummer_unscaled_∇Σ(r, M, a)
end
function Σmean(d::Plummer, r::Real)
    M, a = Mtot(d), scale_radius(d)
    return M / π / (a^2 + r^2) # return 3M / (4π*a^3) * plummer_unscaled_Σmean(r, M, a)
end
function invΣ(d::Plummer, x::Real)
    M, a = Mtot(d), scale_radius(d)
    return a * sqrt( sqrt(M / x / π) / a - 1)
end
function M(d::Plummer, r::Real)
    M, a = Mtot(d), scale_radius(d)
    return a * M * r^3 * sqrt(1 + (r/a)^2) / (a^2 + r^2)^2
end
function ∇M(d::Plummer, r::Real)
    M, a = Mtot(d), scale_radius(d)
    return 3 * a * M * r^2 / (a^2 + r^2)^2 / sqrt(1 + (r/a)^2)
end
# invM is analytic but long and annoying; just use fallback
function Mproj(d::Plummer, r::Real)
    M, a = Mtot(d), scale_radius(d)
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
function Φ(d::Plummer{T}, r::S) where {T, S<:Real}
    U = promote_type(T, S)
    M, a, r = promote(Mtot(d), scale_radius(d), r) 
    return -U(constants.Gvelkpc) * M / sqrt(r^2 + a^2)
end
function ∇Φ(d::Plummer{T}, r::S) where {T, S<:Real}
    U = promote_type(T, S)
    M, a, r = promote(Mtot(d), scale_radius(d), r) 
    return U(constants.Gvelkpc2) * M * r / sqrt((r^2 + a^2)^3)
end
function ∇∇Φ(d::Plummer{T}, r::S) where {T, S<:Real}
    U = promote_type(T, S)
    M, a, r = promote(Mtot(d), scale_radius(d), r) 
    denom = a^2 + r^2
    denom2 = denom^2
    return U(constants.Gvelkpc2) * M * (a^2 - 2r^2) / sqrt(denom2 * denom2 * denom)
end
