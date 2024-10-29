using GalaxyProfiles
using Test

@testset "Plummer" begin
    # Test constructor
    @test Plummer(1.0,1.0) isa Plummer{Float64}
    @test Plummer(1.0,1) isa Plummer{Float64}
    @test Plummer(1.0,1.0f0) isa Plummer{Float64}
    @test Plummer(1.0f0,1.0f0) isa Plummer{Float32}
    types = (Float64, Float32)
    type_labels = ("Float64", "Float32")
    for i in eachindex(types, type_labels)
        @testset "$type_labels[i]" begin
            T = types[i]
            d = Plummer(T(1.0),T(1.0))
            @test @inferred ρ(d,1.0) ≈ 0.04220232731986435
            @test ρ(d,1.0) isa promote_type(T, Float64)
            @test @inferred ρ(d,1.0f0) ≈ 0.04220232731986435f0
            @test ρ(d,1.0f0) isa promote_type(T, Float32)
            @test @inferred invρ(d,0.04220232731986435) ≈ 1.0
            @test invρ(d,0.04220232731986435) isa promote_type(T, Float64)
            @test @inferred invρ(d,0.04220232731986435f0) ≈ 1.0f0
            @test invρ(d,0.04220232731986435f0) isa promote_type(T, Float32)       
            @test @inferred ∇ρ(d,1.0) ≈ -0.10550581829966085
            @test ∇ρ(d,1.0) isa promote_type(T, Float64)
            @test @inferred ∇ρ(d,1.0f0) ≈ -0.10550581829966085f0
            @test ∇ρ(d,1.0f0) isa promote_type(T, Float32)
            @test @inferred ρmean(d,1.0) ≈ 0.0844046546397287
            @test ρmean(d,1.0) isa promote_type(T, Float64)
            @test @inferred ρmean(d,1.0f0) ≈ 0.0844046546397287f0
            @test ρmean(d,1.0f0) isa promote_type(T, Float32)
            @test @inferred invρmean(d,0.0844046546397287) ≈ 1.0
            @test invρmean(d,0.0844046546397287) isa promote_type(T, Float64)
            @test @inferred invρmean(d,0.0844046546397287f0) ≈ 1.0f0
            @test invρmean(d,0.0844046546397287f0) isa promote_type(T, Float32)
            @test @inferred Σ(d,1.0) ≈ 0.07957747154594767
            @test Σ(d,1.0) isa promote_type(T, Float64)
            @test @inferred Σ(d,1.0f0) ≈ 0.07957747154594767f0
            @test Σ(d,1.0f0) isa promote_type(T, Float32)
            @test @inferred ∇Σ(d,1.0) ≈ -0.15915494309189535
            @test ∇Σ(d,1.0) isa promote_type(T, Float64)
            @test @inferred ∇Σ(d,1.0f0) ≈ -0.15915494309189535f0
            @test ∇Σ(d,1.0f0) isa promote_type(T, Float32)
            @test @inferred Σmean(d,1.0) ≈ 0.15915494309189535
            @test Σmean(d,1.0) isa promote_type(T, Float64)
            @test @inferred Σmean(d,1.0f0) ≈ 0.15915494309189535f0
            @test Σmean(d,1.0f0) isa promote_type(T, Float32)
            # ∇Σmean
            @test @inferred invΣ(d,0.07957747154594767) ≈ 1.0
            @test invΣ(d,0.07957747154594767) isa promote_type(T, Float64)
            @test @inferred invΣ(d,0.07957747154594767f0) ≈ 1.0f0
            @test invΣ(d,0.07957747154594767f0) isa promote_type(T, Float32)
            @test @inferred M(d,1.0) ≈ 0.3535533905932738
            @test M(d,1.0) isa promote_type(T, Float64)
            @test @inferred M(d,1.0f0) ≈ 0.3535533905932738f0
            @test M(d,1.0f0) isa promote_type(T, Float32)
            @test @inferred ∇M(d,1.0) ≈ 0.5303300858899106
            @test ∇M(d,1.0) isa promote_type(T, Float64)
            @test @inferred ∇M(d,1.0f0) ≈ 0.5303300858899106f0
            @test ∇M(d,1.0f0) isa promote_type(T, Float32)
            @test @inferred invM(d,0.3535533905932738) ≈ 1.0
            @test invM(d,0.3535533905932738) isa promote_type(T, Float64)
            @test @inferred invM(d,0.3535533905932738f0) ≈ 1.0f0
            @test invM(d,0.3535533905932738f0) isa promote_type(T, Float32)
            @test @inferred M(d,1e9) ≈ Mtot(d) # Test Mtot for Plummer
            @test @inferred Mtot(d) == one(T)  # Test Mtot for Plummer
            @test @inferred Mproj(d,1.0) ≈ 0.5
            @test Mproj(d,1.0) isa promote_type(T, Float64)
            @test @inferred Mproj(d,1.0f0) ≈ 0.5f0
            @test Mproj(d,1.0f0) isa promote_type(T, Float32)
            @test @inferred ∇Mproj(d,1.0) ≈ 0.5
            @test ∇Mproj(d,1.0) isa promote_type(T, Float64)
            @test @inferred ∇Mproj(d,1.0f0) ≈ 0.5f0
            @test ∇Mproj(d,1.0f0) isa promote_type(T, Float32)
            @test @inferred invMproj(d,0.5) ≈ 1.0
            @test invMproj(d,0.5) isa promote_type(T, Float64)
            @test @inferred invMproj(d,0.5f0) ≈ 1.0f0
            @test invMproj(d,0.5f0) isa promote_type(T, Float32)
            @test @inferred Vcirc(d,1.0) ≈ 0.0012331276833655524
            @test Vcirc(d,1.0) isa promote_type(T, Float64)
            @test @inferred Vcirc(d,1.0f0) ≈ 0.0012331276833655524f0
            @test Vcirc(d,1.0f0) isa promote_type(T, Float32)
            @test @inferred Vesc(d,1.0) ≈ 0.0024662553667311048
            @test Vesc(d,1.0) isa promote_type(T, Float64)
            @test @inferred Vesc(d,1.0f0) ≈ 0.0024662553667311048f0
            @test Vesc(d,1.0f0) isa promote_type(T, Float32)
            let d=Plummer(T(10^3), T(10^-2))
                @test @inferred σr(d,0.0,0.0) ≈ 0.267734858583272
                @test σr(d,0.0,0.0) isa promote_type(T, Float64)
                @test @inferred σr(d,0.0f0,0.0f0) ≈ 0.267734858583272f0
                @test σr(d,0.0f0,0.0f0) isa promote_type(T, Float32)
            end
            @test @inferred Φ(d,1.0) ≈ -3.0412077669649876e-6
            @test Φ(d,1.0) isa promote_type(T, Float64)
            @test @inferred Φ(d,1.0f0) ≈ -3.0412077669649876f-6
            @test Φ(d,1.0f0) isa promote_type(T, Float32)
            @test @inferred ∇Φ(d,1.0) ≈ 4.9279415730387377e-23
            @test ∇Φ(d,1.0) isa promote_type(T, Float64)
            @test @inferred ∇Φ(d,1.0f0) ≈ 4.9279415730387377f-23
            @test ∇Φ(d,1.0f0) isa promote_type(T, Float32)
            @test @inferred ∇∇Φ(d,1.0) ≈ -2.4639707865193688e-23
            @test ∇∇Φ(d,1.0) isa promote_type(T, Float64)
            @test @inferred ∇∇Φ(d,1.0f0) ≈ -2.4639707865193688f-23
            @test ∇∇Φ(d,1.0f0) isa promote_type(T, Float32)
            # CDFs and quantiles
            @test @inferred cdf2D(d, 1.0) ≈ 0.5
            @test cdf2D(d, 1.0) isa promote_type(T, Float64)
            @test @inferred cdf2D(d, 1.0f0) ≈ 0.5f0
            @test cdf2D(d, 1.0f0) isa promote_type(T, Float32)
            @test @inferred ccdf2D(d, 1.0) ≈ 0.5
            @test ccdf2D(d, 1.0) isa promote_type(T, Float64)
            @test @inferred ccdf2D(d, 1.0f0) ≈ 0.5f0
            @test ccdf2D(d, 1.0f0) isa promote_type(T, Float32)
            @test @inferred quantile2D(d, 0.5) ≈ 1.0
            @test quantile2D(d, 0.5) isa promote_type(T, Float64)
            @test @inferred quantile2D(d, 0.5f0) ≈ 1.0f0
            @test quantile2D(d, 0.5f0) isa promote_type(T, Float32)
            @test @inferred cquantile2D(d, 0.5) ≈ 1.0
            @test cquantile2D(d, 0.5) isa promote_type(T, Float64)
            @test @inferred cquantile2D(d, 0.5f0) ≈ 1.0f0
            @test cquantile2D(d, 0.5f0) isa promote_type(T, Float32)
            @test @inferred cdf3D(d, 1.0) ≈ 0.3535533905932738
            @test cdf3D(d, 1.0) isa promote_type(T, Float64)
            @test @inferred cdf3D(d, 1.0f0) ≈ 0.3535533905932738f0
            @test cdf3D(d, 1.0f0) isa promote_type(T, Float32)
            @test @inferred ccdf3D(d, 1.0) ≈ 1 - 0.3535533905932738
            @test ccdf3D(d, 1.0) isa promote_type(T, Float64)
            @test @inferred ccdf3D(d, 1.0f0) ≈ 1 - 0.3535533905932738f0
            @test ccdf3D(d, 1.0f0) isa promote_type(T, Float32)
            @test @inferred quantile3D(d, 0.3535533905932738) ≈ 1.0
            @test quantile3D(d, 0.3535533905932738) isa promote_type(T, Float64)
            @test @inferred quantile3D(d, 0.3535533905932738f0) ≈ 1.0f0
            @test quantile3D(d, 0.3535533905932738f0) isa promote_type(T, Float32)
            @test @inferred cquantile3D(d, 1 - 0.3535533905932738) ≈ 1.0
            @test cquantile3D(d, 1 - 0.3535533905932738) isa promote_type(T, Float64)
            @test @inferred cquantile3D(d, 1 - 0.3535533905932738f0) ≈ 1.0f0
            @test cquantile3D(d, 1 - 0.3535533905932738f0) isa promote_type(T, Float32)
            # Other utilities
            @test @inferred GalaxyProfiles.plummer_angular_avalue( T.((-5, 25, 25))... ) ≈
                T(4.324067090092803)
            @test @inferred GalaxyProfiles.plummer_angular_avalue( T.((-5, 25, 27))... ) ≈
                T(1.72144211452034)
            @test @inferred GalaxyProfiles.plummer_a_to_rh( T(20.0) ) ≈ T(26.09532053008214)
            @test GalaxyProfiles.plummer_a_to_rh( T(20.0) ) isa T
            @test @inferred GalaxyProfiles.plummer_rh_to_a( T(26.09532053008214) ) ≈ T(20.0)
            @test GalaxyProfiles.plummer_rh_to_a( T(26.09532053008214) ) isa T
        end
    end
end


# @testset "Float64" begin
#     @test Plummer(1.0,1.0) isa Plummer{Float64}
#     @test Plummer(1.0,1) isa Plummer{Float64}
#     @test Plummer(1.0,1.0f0) isa Plummer{Float64}
#     d = Plummer(1.0,1.0)
#     @test @inferred ρ(d,1.0) ≈ 0.04220232731986435
#     @test ρ(d,1.0) isa Float64
#     @test @inferred ρ(d,1.0f0) ≈ 0.04220232731986435f0
#     @test ρ(d,1.0f0) isa Float64
#     @test @inferred invρ(d,0.04220232731986435) ≈ 1.0
#     @test invρ(d,0.04220232731986435) isa Float64
#     @test @inferred invρ(d,0.04220232731986435) ≈ 1.0f0
#     @test invρ(d,0.25f0) isa Float64        
#     @test @inferred ∇ρ(d,1.0) ≈ -0.10550581829966085
#     @test ∇ρ(d,1.0) isa Float64
#     @test @inferred ∇ρ(d,1.0f0) ≈ -0.10550581829966085f0
#     @test ∇ρ(d,1.0f0) isa Float64
#     @test @inferred ρmean(d,1.0) ≈ 0.0844046546397287
#     @test ρmean(d,1.0) isa Float64
#     @test @inferred ρmean(d,1.0f0) ≈ 0.0844046546397287f0
#     @test ρmean(d,1.0f0) isa Float64
#     @test @inferred invρmean(d,0.0844046546397287) ≈ 1.0
#     @test invρmean(d,0.0844046546397287) isa Float64
#     @test @inferred invρmean(d,0.0844046546397287f0) ≈ 1.0f0
#     @test invρmean(d,0.0844046546397287f0) isa Float64
#     @test @inferred Σ(d,1.0) ≈ 0.07957747154594767
#     @test Σ(d,1.0) isa Float64
#     @test @inferred Σ(d,1.0f0) ≈ 0.07957747154594767f0
#     @test Σ(d,1.0f0) isa Float64
#     @test @inferred ∇Σ(d,1.0) ≈ -0.15915494309189535
#     @test ∇Σ(d,1.0) isa Float64
#     @test @inferred ∇Σ(d,1.0f0) ≈ -0.15915494309189535f0
#     @test ∇Σ(d,1.0f0) isa Float64
#     @test @inferred Σmean(d,1.0) ≈ 0.15915494309189535
#     @test Σmean(d,1.0) isa Float64
#     @test @inferred Σmean(d,1.0f0) ≈ 0.15915494309189535f0
#     @test Σmean(d,1.0f0) isa Float64
#     # ∇Σmean
#     @test @inferred invΣ(d,0.07957747154594767) ≈ 1.0
#     @test invΣ(d,0.07957747154594767) isa Float64
#     @test @inferred invΣ(d,0.07957747154594767f0) ≈ 1.0f0
#     @test invΣ(d,0.07957747154594767f0) isa Float64
#     @test @inferred M(d,1.0) ≈ 0.3535533905932738
#     @test M(d,1.0) isa Float64
#     @test @inferred M(d,1.0f0) ≈ 0.3535533905932738f0
#     @test M(d,1.0f0) isa Float64
#     @test @inferred ∇M(d,1.0) ≈ 0.5303300858899106
#     @test ∇M(d,1.0) isa Float64
#     @test @inferred ∇M(d,1.0f0) ≈ 0.5303300858899106f0
#     @test ∇M(d,1.0f0) isa Float64
#     @test @inferred invM(d,0.3535533905932738) ≈ 1.0
#     @test invM(d,0.3535533905932738) isa Float64
#     @test @inferred invM(d,0.3535533905932738f0) ≈ 1.0f0
#     @test invM(d,0.3535533905932738f0) isa Float64
#     @test @inferred M(d,1e9) ≈ Mtot(d) # Test Mtot for Plummer
#     @test @inferred Mtot(d) == 1.0     # Test Mtot for Plummer
#     @test @inferred Mproj(d,1.0) ≈ 0.5
#     @test Mproj(d,1.0) isa Float64
#     @test @inferred Mproj(d,1.0f0) ≈ 0.5f0
#     @test Mproj(d,1.0f0) isa Float64
#     @test @inferred ∇Mproj(d,1.0) ≈ 0.5
#     @test ∇Mproj(d,1.0) isa Float64
#     @test @inferred ∇Mproj(d,1.0f0) ≈ 0.5f0
#     @test ∇Mproj(d,1.0f0) isa Float64
#     @test @inferred invMproj(d,0.5) ≈ 1.0
#     @test invMproj(d,0.5) isa Float64
#     @test @inferred invMproj(d,0.5f0) ≈ 1.0f0
#     @test invMproj(d,0.5f0) isa Float64
#     @test @inferred Vcirc(d,1.0) ≈ 0.0012331276833655524
#     @test Vcirc(d,1.0) isa Float64
#     @test @inferred Vcirc(d,1.0f0) ≈ 0.0012331276833655524f0
#     @test Vcirc(d,1.0f0) isa Float64
#     @test @inferred Vesc(d,1.0) ≈ 0.0024662553667311048
#     @test Vesc(d,1.0) isa Float64
#     @test @inferred Vesc(d,1.0f0) ≈ 0.0024662553667311048f0
#     @test Vesc(d,1.0f0) isa Float64
#     @test @inferred Φ(d,1.0) ≈ -3.0412077669649876e-6
#     @test Φ(d,1.0) isa Float64
#     @test @inferred Φ(d,1.0f0) ≈ -3.0412077669649876f-6
#     @test Φ(d,1.0f0) isa Float64
#     @test @inferred ∇Φ(d,1.0) ≈ 4.9279415730387377e-23
#     @test ∇Φ(d,1.0) isa Float64
#     @test @inferred ∇Φ(d,1.0f0) ≈ 4.9279415730387377f-23
#     @test ∇Φ(d,1.0f0) isa Float64
#     @test @inferred ∇∇Φ(d,1.0) ≈ -2.4639707865193688e-23
#     @test ∇∇Φ(d,1.0) isa Float64
#     @test @inferred ∇∇Φ(d,1.0f0) ≈ -2.4639707865193688f-23
#     @test ∇∇Φ(d,1.0f0) isa Float64
#     # CDFs and quantiles
#     @test @inferred cdf2D(d, 1.0) ≈ 0.5
#     @test cdf2D(d, 1.0) isa Float64
#     @test @inferred cdf2D(d, 1.0f0) ≈ 0.5f0
#     @test cdf2D(d, 1.0f0) isa Float64
#     @test @inferred ccdf2D(d, 1.0) ≈ 0.5
#     @test ccdf2D(d, 1.0) isa Float64
#     @test @inferred ccdf2D(d, 1.0f0) ≈ 0.5f0
#     @test ccdf2D(d, 1.0f0) isa Float64
#     @test @inferred quantile2D(d, 0.5) ≈ 1.0
#     @test quantile2D(d, 0.5) isa Float64
#     @test @inferred quantile2D(d, 0.5f0) ≈ 1.0f0
#     @test quantile2D(d, 0.5f0) isa Float64
#     @test @inferred cquantile2D(d, 0.5) ≈ 1.0
#     @test cquantile2D(d, 0.5) isa Float64
#     @test @inferred cquantile2D(d, 0.5f0) ≈ 1.0f0
#     @test cquantile2D(d, 0.5f0) isa Float64
#     @test @inferred cdf3D(d, 1.0) ≈ 0.3535533905932738
#     @test cdf3D(d, 1.0) isa Float64
#     @test @inferred cdf3D(d, 1.0f0) ≈ 0.3535533905932738f0
#     @test cdf3D(d, 1.0f0) isa Float64
#     @test @inferred ccdf3D(d, 1.0) ≈ 1 - 0.3535533905932738
#     @test ccdf3D(d, 1.0) isa Float64
#     @test @inferred ccdf3D(d, 1.0f0) ≈ 1 - 0.3535533905932738f0
#     @test ccdf3D(d, 1.0f0) isa Float64
#     @test @inferred quantile3D(d, 0.3535533905932738) ≈ 1.0
#     @test quantile3D(d, 0.3535533905932738) isa Float64
#     @test @inferred quantile3D(d, 0.3535533905932738f0) ≈ 1.0f0
#     @test quantile3D(d, 0.3535533905932738f0) isa Float64
#     @test @inferred cquantile3D(d, 1 - 0.3535533905932738) ≈ 1.0
#     @test cquantile3D(d, 1 - 0.3535533905932738) isa Float64
#     @test @inferred cquantile3D(d, 1 - 0.3535533905932738f0) ≈ 1.0f0
#     @test cquantile3D(d, 1 - 0.3535533905932738f0) isa Float64
# end
