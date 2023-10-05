# Overall need to decide on an API for mass profiles vs light profiles. Just dumping some code here for now.
# I think I should favor mass/light style implementations and then provide some sort of conversion utility for surface brightness profiles. 

# struct Sersic{T<:Real} # <: AbstractSurfaceDensity
#     μ_e::T
#     r_e::T
#     n::T
#     q::T
#     # Sersic{T}(μ_e::T,r_e::T,n::T,q::T) where {T} = new{T}(μ_e,r_e,n,q)
# end
# Sersic(μ_e::Real, r_e::Real, n::Real, q::Real) = Sersic(promote(μ_e,r_e,n,q)...)

struct Sersic{T<:Real} # <: AbstractSurfaceDensity
    Σ0::T
    r_e::T
    n::T
    q::T
    # Sersic{T}(μ_e::T,r_e::T,n::T,q::T) where {T} = new{T}(μ_e,r_e,n,q)
end
Sersic(Σ0::Real, r_e::Real, n::Real, q::Real) = Sersic(promote(Σ0,r_e,n,q)...)

#### Parameters
scale_radius(d::Sersic) = d.r_e
params(d::Sersic) = (d.Σ0,d.r_e,d.n,d.q)

b_n(n::T) where T<:Real = gamma_inc_inv(2n, T(1//2), T(1//2))
b_n(d::Sersic) = b_n(d.n)
f_n(n::Real) = (bn = b_n(n); gamma(2n) * n * exp(bn) / bn^(2n)) # this f_n from graham 2005, equation 8
f_n(d::Sersic) = f_n(d.n) # this f_n from graham 2005, equation 8


#### Conversions
Base.convert(::Type{Sersic{T}}, d::Sersic) where T = Sersic(convert(T,d.μ_e), convert(T,d.r_e), convert(T,d.n), convert(T,d.q))
Base.convert(::Type{Sersic{T}}, d::Sersic{T}) where {T<:Real} = d

#### Evaluation
function Σ(d::Sersic, r::Real) 
    Σ0, r_e, n, q = params(d)
    Σ0 * exp( -b_n(d) * ( (r/r_e)^inv(n) - 1))
end
function ∇Σ(d::Sersic, r::Real)
    Σ0, r_e, n, q = params(d)
    -Σ(d, r) * b_n(d) * (r/r_e)^(inv(n) - 1) / n / r_e
end
# function Σmean(d::Sersic, r::Real)
#     Σ0, r_e, n, q = params(d)
#     bn = b_n(d)
#     rfac = bn * (r/r_e)^inv(n)
#     Σ0 * exp(bn) * n * r^2 * gamma(2n, rfac) * (rfac)^(-2n)
# end
function invΣ(d::Sersic, x::Real)
    Σ0, r_e, n, q = params(d)
    x<=0 ? throw(DomainError(x, "x must be greater 0")) : r_e * (-log(x/Σ0) / b_n(d) + 1)^n
end
function Mproj(d::Sersic{T}, r::S) where {T, S<:Real}
    U = promote_type(T, S)
    isinf(r) && (return U(Mtot(d)))
    Σ0, r_e, n, q = params(d)
    bn = b_n(d)
    # rfac = bn * (r/r_e)^inv(n)
    # U(2π) * Σ0 * n * r_e^2 * exp(bn) * bn^-2n * (gamma(2n) - gamma(2n, bn*(r/r_e)^inv(n)))
    U(2π) * Σ0 * n * r_e^2 * exp(bn) * bn^-2n * gamma_inc(2n, bn*(r/r_e)^inv(n))[1] * gamma(2n)
end
function ∇Mproj(d::Sersic{T}, r::S) where {T, S<:Real}
    U = promote_type(T, S)
    Σ0, r_e, n, q = params(d)
    U(2π) * r * Σ(d, r)
end
# invMproj
function Mtot(d::Sersic{T}) where T
    Σ0, r_e, n, q = params(d)
    bn = b_n(d)
    T(2π) * Σ0 * n * r_e^2 * exp(bn) * bn^-2n * gamma(2n)
end
