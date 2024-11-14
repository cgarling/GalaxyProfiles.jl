# Test common fallback functions by defining a new struct that matches
# a profile already defined and subset of necessary methods on that struct,
# and testing the common fallbacks against the profile-specific methods;
# both should give the same answer to high precision.

using GalaxyProfiles
import GalaxyProfiles: ρ, params, scale_radius, Mtot, ∇M, AbstractSurfaceDensity, Σ # Explicitly import to extend methods
using Test

# Define new Plummer struct -- implement only subset of methods
# to test the fallbacks against the methods specialized on the Plummer type
struct Plummer2{T <: Real} <: GalaxyProfiles.AbstractDensity{T}
    M::T
    a::T
end
Plummer2(M::Real, a::Real) = Plummer2(promote(M,a)...)
params(d::Plummer2) = (d.M, d.a)
scale_radius(d::Plummer2) = d.a
Mtot(d::Plummer2) = d.M
ρ(d::Plummer2, r::Real) = ρ(Plummer(params(d)...), r)
∇M(d::Plummer2, r::Real) = ∇M(Plummer(params(d)...), r) # Required to test ∇∇Φ

@testset "Common Methods -- Densities" begin
    # Construct instances and run tests
    types = (Float64, Float32)
    type_labels = ("Float64", "Float32")
    for i in eachindex(types, type_labels)
        tl = type_labels[i]
        @testset "$tl" begin
            T = types[i]
            mass, a = T(1e3), T(2)
            r = T(1//2)
            p = Plummer(mass, a)
            p2 = Plummer2(mass, a)

            # Test common methods
            @test ρ(p, r) ≈ ρ(p2, r)
            @test invρ(p, ρ(p,r)) ≈ invρ(p2, ρ(p2,r))
            # @test ∇ρ(p, r) ≈ ∇ρ(p2, r)
            @test ρmean(p, r) ≈ ρmean(p2, r)
            @test invρmean(p, ρmean(p,r)) ≈ invρmean(p2, ρmean(p2,r))
            @test Σ(p, r) ≈ Σ(p2, r)
            # @test ∇Σ(p, r) ≈ ∇Σ(p2, r)
            @test Σmean(p, r) ≈ Σmean(p2, r)
            @test invΣ(p, Σ(p,r)) ≈ invΣ(p2, Σ(p2,r))
            @test M(p, r) ≈ M(p2, r) # Plummer uses the fallback here so not a real test perhaps...
            # @test ∇M(p, r) ≈ ∇M(p2, r)
            @test invM(p, M(p,r)) ≈ invM(p2, M(p2,r))
            # @test Mtot(p, r) ≈ Mtot(p2, r)
            @test Mproj(p, r) ≈ Mproj(p2, r)
            # @test ∇Mproj(p, r) ≈ ∇Mproj(p2, r)
            @test invMproj(p, Mproj(p,r)) ≈ invMproj(p2, Mproj(p2,r))
            # dynamical_time always uses the common function
            @test cdf2D(p, r) ≈ cdf2D(p2, r)
            # ccdf2D always uses common function
            @test cdf3D(p, r) ≈ cdf3D(p2, r)
            # ccdf3D always uses common function
            # @test quantile2D(p, r) ≈ quantile2D(p2, r)
            # @test cquantile2D(p, r) ≈ quantile2D(p2, r)
            # @test quantile3D(p, r) ≈ quantile3D(p2, r)
            # @test cquantile3D(p, r) ≈ quantile3D(p2, r)
            # Vcirc uses common function
            # Vesc uses common function
            # No Vmax common function
            @test σr(p, r, zero(T)) ≈ σr(p2, r, zero(T))
            @test σlos(p, r, zero(T)) ≈ σlos(p2, r, zero(T))
            @test Φ(p, r) ≈ Φ(p2, r)
            @test ∇Φ(p, r) ≈ ∇Φ(p2, r)
            @test ∇∇Φ(p, r) ≈ ∇∇Φ(p2, r)
        end
    end
end


# Define new ExponentialDisk struct -- implement only subset of methods
# to test the fallbacks against the methods specialized on the ExponentialDisk type
struct ExponentialDisk2{T <: Real} <: AbstractSurfaceDensity{T}
    Σ0::T
    rs::T
end
ExponentialDisk2(Σ0::Real, rs::Real) = ExponentialDisk2(promote(Σ0,rs)...)
params(d::ExponentialDisk2) = (d.Σ0, d.rs)
scale_radius(d::ExponentialDisk2) = d.rs
Mtot(d::ExponentialDisk2) = Mtot(ExponentialDisk(params(d)...))
Σ(d::ExponentialDisk2, r::Real) = Σ(ExponentialDisk(params(d)...), r)

@testset "Common Methods -- Surface Densities" begin
    # Construct instances and run tests
    types = (Float64, Float32)
    type_labels = ("Float64", "Float32")
    for i in eachindex(types, type_labels)
        tl = type_labels[i]
        @testset "$tl" begin
            T = types[i]
            Σ0, rs = T(1e2), T(2)
            r = T(1//2)
            p = ExponentialDisk(Σ0, rs)
            p2 = ExponentialDisk2(Σ0, rs)

            # Test common methods
            # @test ∇Σ(p, r) ≈ ∇Σ(p2, r)
            @test Σmean(p, r) ≈ Σmean(p2, r)
            @test invΣ(p, Σ(p,r)) ≈ invΣ(p2, Σ(p2,r))
            @test Mproj(p, r) ≈ Mproj(p2, r)
            # @test ∇Mproj(p, r) ≈ ∇Mproj(p2, r)
            @test invMproj(p, Mproj(p,r)) ≈ invMproj(p2, Mproj(p2,r))
            @test cdf2D(p, r) ≈ cdf2D(p2, r)
        end
    end
end
