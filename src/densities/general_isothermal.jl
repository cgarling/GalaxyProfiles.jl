
"""
    GeneralIsothermal(ρ0::Real, rs::Real, α::Real)
    GeneralIsothermal(ρ0::Unitful.Density, rs::Unitful.Length, α::Real)
    GeneralIsothermal(rs::Real, α::Real, M::Real, Rmax::Real)
    GeneralIsothermal(rs::Unitful.Length, α::Real, M::Unitful.Mass, Rmax::Unitful.Length)

Type describing general isothermal density profiles with scale radius `rs`, power-law index `α`, and density at `rs` of `ρ0`. The density profile is

```math
\\rho(r) = \\rho_0 \\times \\left( \\frac{r}{R_s} \\right)^{-\\alpha}
```

The fields of `GeneralIsothermal` are `ρ0, rs, α`. The default units of `GeneralIsothermal` are `[ρ0] = [Msun/kpc^3], [r, rs] = [kpc], [M] = [Msun]` when you construct `GeneralIsothermal` with `Real`s, like `Float64`. This is important for quantities like [`Vcirc`](@ref), [`Vesc`](@ref), [`Φ`](@ref), [`∇Φ`](@ref), and [`∇∇Φ`](@ref) which involve `G`; these will give incorrect results if the fields of `GeneralIsothermal` or the provided `r` are in different units. If you construct `GeneralIsothermal` with `Unitful` quantities, they will be internally converted. 

Since the total mass of the `GeneralIsothermal` profile is undefined when ``\\alpha \\geq 3``, we define the the potential [`Φ`](@ref) to be 0 at `rs` for *all* instances of `GeneralIsothermal`.

The following methods are specialized on this type:
 - [`ρ`](@ref), [`invρ`](@ref), [`∇ρ`](@ref), [`Σ`](@ref), [`∇Σ`](@ref), [`invΣ`](@ref), [`M`](@ref), [`∇M`](@ref), [`invM`](@ref), [`Mproj`](@ref), [`∇Mproj`](@ref), [`invMproj`](@ref), [`Vcirc`](@ref), [`Vesc`](@ref), [`Φ`](@ref), [`∇Φ`](@ref), [`∇∇Φ`](@ref)

# See also
 - [`SIS`](@ref)
"""
struct GeneralIsothermal{T <: Real} <: AbstractSurfaceDensity
    ρ0::T
    rs::T
    α::T
end
GeneralIsothermal(ρ0::Real, rs::Real, α::Real) = GeneralIsothermal(promote(ρ0,rs,α)...)
generalisothermal_from_M_Rmax(rs::T, α::T, M::T, Rmax::T) where {T<:Real} =
    α > 3 ?
    throw(DomainError(α,"Enclosed mass is only finite for α<3.")) :
    GeneralIsothermal(M * (3 - α) / (rs^α * Rmax^(3-α) * 4 * π), rs, α)
GeneralIsothermal(rs::Real, α::Real, M::Real, Rmax::Real) = generalisothermal_from_M_Rmax(promote(rs,α,M,Rmax)...)

"""
    SIS(ρ0::Real, rs::Real)
    SIS(ρ0::Unitful.Density, rs::Unitful.Length)
    SIS(rs::Real, M::Real, Rmax::Real)
    SIS(rs::Unitful.Length, M::Unitful.Mass, Rmax::Unitful.Length)

Convenience function to construct a singular isothermal sphere; i.e., a [`GeneralIsothermal`](@ref) with `α=2`.
"""
SIS(ρ0::Real, rs::Real) = GeneralIsothermal(ρ0, rs, 2)
SIS(rs::Real, M::Real, Rmax::Real) = GeneralIsothermal(rs, 2, M, Rmax)

#### Parameters

params(d::GeneralIsothermal) = (d.ρ0, d.rs, d.α)
scale_radius(d::GeneralIsothermal) = d.rs

#### Evaluation

function ρ(d::GeneralIsothermal, r::Real)
    ρ0, rs, α = params(d)
    return ρ0 * (r/rs)^-α
end
function invρ(d::GeneralIsothermal, x::Real)
    ρ0, rs, α = params(d)
    return rs * (x/ρ0)^(-1/α)
end
function ∇ρ(d::GeneralIsothermal,r::Real)
    ρ0, rs, α = params(d)
    return -α * ρ0 * (r/rs)^-(α+1) / rs
end
# ρmean
# invρmean
function Σ(d::GeneralIsothermal{T}, r::S) where {T, S <: Real}
    # This is the abel integral. Because it integrates to infinite along the line of sight, the integral of Σ,
    # e.g. quadgk(x->2π*x*Σ(d,x),0,R) will not equal 10 at Rmax, if you used the total mass enclosed within
    # See Equation 2.59 in Binney Galactic Dynamics 2E, page 81 and appendix c.2 on page 798 about the factorials

    U = promote_type(T, S)
    ρ0, rs, α = params(d)
    ρ0, rs, α, r = promote(ρ0, rs, α, r)
    α>=3 && throw(DomainError(α, "Σ is only finite for α<3."))
    return sqrt(U(π)) * ρ0 * rs^α * gamma((α-1)/2) / gamma(α/2) / r^(α-1)
end
function ∇Σ(d::GeneralIsothermal{T}, r::S) where {T, S <: Real}
    U = promote_type(T, S)
    ρ0, rs, α = params(d)
    ρ0, rs, α, r = promote(ρ0, rs, α, r)
    α>=3 && throw(DomainError(α, "Σ is only finite for α<3."))
    return -sqrt(U(π)) * ρ0 * rs^α * (α-1) * gamma((α-1)/2) / gamma(α/2) / r^α
end
# Σmean
function invΣ(d::GeneralIsothermal{T}, x::S) where {T, S <: Real}
    U = promote_type(T, S)
    ρ0, rs, α = params(d)
    ρ0, rs, α, x = promote(ρ0, rs, α, x)
    if α>=3
        throw(DomainError(α, "Σ is only finite for α<3."))
    elseif x<=0
        throw(DomainError(x, "invΣ only valid for x>0."))
    else
        return (sqrt(U(π)) * ρ0 * rs^α * gamma((α-1)/2) / gamma(α/2) / x)^(α-1)
    end
end
function M(d::GeneralIsothermal{T}, r::S) where {T, S <: Real}
    ρ0, rs, α = params(d)
    α >= 3  && throw(DomainError(α, "Enclosed mass is only finite for α<3."))
    return ρ0 * rs^α * r^(3-α) * 4 * π / (3 - α)
end
function ∇M(d::GeneralIsothermal{T}, r::S) where {T, S <: Real}
    ρ0, rs, α = params(d)
    α >= 3 && throw(DomainError(α,"Enclosed mass is only finite for α<3."))
    return ρ0 * rs^α * r^(2-α) * 4 * π
end
function invM(d::GeneralIsothermal{T}, x::S) where {T, S <: Real}
    ρ0,rs,α = params(d)
    ρ0, rs, α, x = promote(ρ0, rs, α, x)
    α >= 3 && throw(DomainError(α, "Enclosed mass is only finite for α<3."))
    return (x * (3 - α) / (ρ0 * rs^α * 4 * π))^(1/(3-α))
end
# Mtot
function Mproj(d::GeneralIsothermal{T}, r::S) where {T, S <: Real}
    U = promote_type(T, S)
    ρ0, rs, α = params(d)
    ρ0, rs, α, r = promote(ρ0, rs, α, r)
    α >= 3 && throw(DomainError(α, "Mproj is only finite for α<3."))
    return 2 * U(π^(3/2)) * r^(3-α) * rs^α * ρ0 * gamma((α-1)/2) / (3-α) / gamma(α/2)
end
function ∇Mproj(d::GeneralIsothermal{T}, r::S) where {T, S <: Real}
    U = promote_type(T,S)
    ρ0,rs,α = params(d)
    ρ0, rs, α, r = promote(ρ0, rs, α, r)
    α >= 3 && throw(DomainError(α, "Σ is only finite for α<3."))
    return 2 * U(π^(3/2)) * r^(2-α) * rs^α * ρ0 * gamma((α-1)/2) / gamma(α/2)
end
function invMproj(d::GeneralIsothermal{T}, x::S) where {T, S <: Real}
    U = promote_type(T, S)
    ρ0, rs, α = params(d)
    ρ0, rs, α, x = promote(ρ0, rs, α, x)
    if α >= 3
        throw(DomainError(α, "Σ is only finite for α<3."))
    elseif x<=0
        throw(DomainError(x, "invΣ only valid for x>0."))
    else
        (2 * U(π^(3/2)) * rs^α * ρ0 * gamma((α-1)/2) / (3-α) / gamma(α/2) / x)^(3-α)
    end
end
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
# cquantile
# Equation 2.61 on page 81 of Binney Galactic Dynamics 2E
function Vcirc(d::GeneralIsothermal{T}, r::S) where {T, S <: Real}
    U = promote_type(T, S)
    ρ0, rs, α = params(d)
    ρ0, rs, α, r = promote(ρ0, rs, α, r)
    # α >= 3 ? throw(DomainError(α,"Circular velocity is only finite for α<3.")) : sqrt( 4π * constants.Gvelkpc * ρ0 * rs^α * r^(2-α) / (3 - α))
    α >= 3 && throw(DomainError(α, "Circular velocity is only finite for α<3."))
    if α == 2
        sqrt(4 * U(π) * U(constants.Gvelkpc) * ρ0 * rs^α)
    else
        sqrt(4 * U(π) * U(constants.Gvelkpc) * ρ0 * rs^α * r^(2-α) / (3 - α))
    end
end
# Equation 2.63 on page 81 of Binney Galactic Dynamics 2E
function Vesc(d::GeneralIsothermal{T} ,r::S) where {T, S <: Real}
    U = promote_type(T, S)
    ρ0, rs, α = params(d)
    ρ0, rs, α, r = promote(ρ0, rs, α, r)
    # α <= 2 && throw(DomainError(α,"Escape velocity is only finite for α>2."))
    # 8 * U(π) * U(constants.Gvelkpc) * ρ0 * rs^α * r^(2-α) / (3 - α) / (α - 2)
    if 2 < α < 3 
        8 * U(π) * U(constants.Gvelkpc) * ρ0 * rs^α * r^(2-α) / (3 - α) / (α - 2)
    else
        throw(DomainError(α, "Escape velocity is only finite for α>2."))
    end
end
# Equation 2.62 on page 81 of Binney Galactic Dynamics 2E
function Φ(d::GeneralIsothermal{T}, r::S) where {T, S <: Real}
    U = promote_type(T, S)
    ρ0, rs, α = params(d)
    ρ0, rs, α, r = promote(ρ0, rs, α, r)
    # this is correct but slow
    # α == 2 ? Vcirc(d,r)^2 * log(r/rs) : (Vcirc(d,rs)^2 - Vcirc(d,r)^2)/(α-2)
    α >= 3 && throw(DomainError(α, "Potential is only finite for α<3."))
    if α == 2
        return 4 * U(π) * U(constants.Gvelkpc) * ρ0 * rs^α * log(r/rs)# / (3 - α) 
    else
        vrs = 4 * U(π) * U(constants.Gvelkpc) * ρ0 * rs^2 / (3 - α)
        vr = 4 * U(π) * U(constants.Gvelkpc) * ρ0 * rs^α * r^(2-α) / (3 - α)
        return (vrs-vr)/(α-2)
        # 4π * constants.Gvelkpc * ρ0 * (rs^2-(rs^α * r^(2-α))) / ((3 - α)*(α-2))
    end
end
function ∇Φ(d::GeneralIsothermal{T}, r::S) where {T, S <: Real} # returns in km/s^2
    # (PhysicalConstants.CODATA2018.G |> ua.kpc^2*u.km/u.s^2/ua.Msun) * (ua.Msun/ua.kpc^3) * ua.kpc^2 / ua.kpc = u.km/u.s^2
    # (PhysicalConstants.CODATA2018.G |> ua.kpc*u.km^2/u.s^2/ua.Msun) * (ua.Msun/ua.kpc^3) * ua.kpc^2 / ua.kpc = u.km^2/u.s^2/ua.kpc
    U = promote_type(T, S)
    ρ0, rs, α = params(d)
    ρ0, rs, α, r = promote(ρ0, rs, α, r)
    # this is correct but slow
    # α == 2 ? Vcirc(d,r)^2 * log(r/rs) : (Vcirc(d,rs)^2 - Vcirc(d,r)^2)/(α-2)
    α >= 3 && throw(DomainError(α, "Potential is only finite for α<3."))
    if α == 2
        # 4π * constants.Gvelkpc2 * ρ0 * rs^2 / r# / (3 - α) 
        return 4 * U(π) * U(constants.Gvelkpc2) * ρ0 * rs^2 / r # / (3 - α) 
    else
        # -4π * constants.Gvelkpc2 * ρ0 * r^(1-α) * rs^α * (2-α) / ((3-α)*(α-2))
        return -4 * U(π) * U(constants.Gvelkpc2) * ρ0 * r^(1-α) * rs^α * (2-α) / ((3-α)*(α-2))
    end
end
function ∇∇Φ(d::GeneralIsothermal{T}, r::S) where {T, S <: Real} # returns in km/s^2/kpc
    U = promote_type(T, S)
    ρ0, rs, α = params(d)
    ρ0, rs, α, r = promote(ρ0, rs, α, r)
    α >= 3 && throw(DomainError(α, "Potential is only finite for α<3."))
    if α == 2
        -4 * U(π) * U(constants.Gvelkpc2) * ρ0 * rs^2 / r^2 # / (3 - α) 
    else
        -4 * U(π) * U(constants.Gvelkpc2) * ρ0 * rs^α * (2-α) * (1-α) / ((3-α)*(α-2)) / r^α
    end
end
