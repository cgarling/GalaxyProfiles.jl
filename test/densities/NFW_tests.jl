# from galpy.potential import NFWPotential, tdyn
# from astropy import units as u
# # Galpy uses NFW amplitude in units of mass rather than density, so need to correct for that
# p = NFWPotential(amp=1e3 * u.Msun / u.kpc**3 * (4 * np.pi * (10*u.kpc)**3), a=10 * u.kpc)
# print(tdyn(p, 10*u.kpc))
# # 60.130773 Gyr
# print(tdyn(p, 2*10*u.kpc))
# # 113.72905 Gyr

using GalaxyProfiles
using Test

@testset "NFW" begin
    @testset "Float64" begin
        @test NFW(1.0,1.0) isa NFW{Float64}
        @test NFW(1.0,1) isa NFW{Float64}
        @test NFW(1.0,1.0f0) isa NFW{Float64}
        d = NFW(1.0,1.0)
        @test @inferred ρ(d,1.0) ≈ 0.25
        @test ρ(d,1.0) isa Float64
        @test @inferred ρ(d,1.0f0) ≈ 0.25f0
        @test ρ(d,1.0f0) isa Float64
        @test @inferred invρ(d,0.25) ≈ 1.0
        @test invρ(d,0.25) isa Float64
        @test @inferred invρ(d,0.25f0) ≈ 1.0f0
        @test invρ(d,0.25f0) isa Float64        
        @test @inferred ∇ρ(d,1.0) ≈ -0.375
        @test ∇ρ(d,1.0) isa Float64
        @test @inferred ∇ρ(d,1.0f0) ≈ -0.375f0
        @test ∇ρ(d,1.0f0) isa Float64
        @test @inferred ρmean(d,1.0) ≈ 0.5794415416798359
        @test ρmean(d,1.0) isa Float64
        @test @inferred ρmean(d,1.0f0) ≈ 0.5794415416798359f0
        @test ρmean(d,1.0f0) isa Float64
        @test @inferred invρmean(d,0.2) ≈ 1.795108201392415
        @test invρmean(d,0.2) isa Float64
        @test @inferred invρmean(d,0.2f0) ≈ 1.795108201392415f0
        @test invρmean(d,0.2f0) isa Float64
        @test @inferred Σ(d,1.0) ≈ 2//3
        @test Σ(d,1.0) isa Float64
        @test @inferred Σ(d,1.0f0) ≈ Float32(2//3)
        @test Σ(d,1.0f0) isa Float64
        # ∇Σ
        @test @inferred Σmean(d,1.0) ≈ 1.2274112777602189
        @test Σmean(d,1.0) isa Float64
        @test @inferred Σmean(d,1.0f0) ≈ 1.2274112777602189f0
        @test Σmean(d,1.0f0) isa Float64
        # ∇Σmean
        @test @inferred invΣ(d,0.5) ≈ 1.259435551338637
        @test invΣ(d,0.5) isa Float64
        @test @inferred invΣ(d,0.5f0) ≈ 1.259435551338637f0
        @test invΣ(d,0.5f0) isa Float64
        @test @inferred M(d,1.0) ≈ 2.4271590540348216
        @test M(d,1.0) isa Float64
        @test @inferred M(d,1.0f0) ≈ 2.4271590540348216f0
        @test M(d,1.0f0) isa Float64
        @test @inferred ∇M(d,1.0) ≈ π
        @test ∇M(d,1.0) isa Float64
        @test @inferred ∇M(d,1.0f0) ≈ Float32(π)
        @test ∇M(d,1.0f0) isa Float64
        @test @inferred invM(d,2.4271590540348216) ≈ 1.0
        @test invM(d,2.4271590540348216) isa Float64
        @test @inferred invM(d,2.4271590540348216f0) ≈ 1.0f0
        @test invM(d,2.4271590540348216f0) isa Float64
        # Mtot
        @test @inferred Mproj(d,1.0) ≈ 3.8560262531447647
        @test Mproj(d,1.0) isa Float64
        @test @inferred Mproj(d,1.0f0) ≈ 3.8560262531447647f0
        @test Mproj(d,1.0f0) isa Float64
        @test @inferred ∇Mproj(d,1.0) ≈ 4.1887902047863905
        @test ∇Mproj(d,1.0) isa Float64
        @test @inferred ∇Mproj(d,1.0f0) ≈ 4.1887902047863905f0
        @test ∇Mproj(d,1.0f0) isa Float64
        @test @inferred invMproj(d,0.2) ≈ 0.11585086460387914
        @test invMproj(d,0.2) isa Float64
        @test @inferred invMproj(d,0.2f0) ≈ 0.11585086460387914f0
        @test invMproj(d,0.2f0) isa Float64
        @test @inferred dynamical_time(d,1.0) ≈ 4.753754967831782e11
        @test dynamical_time(d,1.0) isa Float64
        @test @inferred dynamical_time(d,1.0f0) ≈ 4.753754967831782f11
        @test dynamical_time(d,1.0f0) isa Float64
        @test @inferred Vcirc(d,1.0) ≈ 0.003230945727279133
        @test Vcirc(d,1.0) isa Float64
        @test @inferred Vcirc(d,1.0f0) ≈ 0.003230945727279133f0
        @test Vcirc(d,1.0f0) isa Float64
        @test @inferred Vesc(d,1.0) ≈ 0.008655919418653363
        @test Vesc(d,1.0) isa Float64
        @test @inferred Vesc(d,1.0f0) ≈ 0.008655919418653363f0
        @test Vesc(d,1.0f0) isa Float64
        @test @inferred all(Vmax(d) .≈ (0.003418455956363109, 2.1625815870646097))
        @test Vmax(d) isa NTuple{2,Float64}
        @test @inferred Φ(d,1.0) ≈ -3.7462470491110184e-5
        @test Φ(d,1.0) isa Float64
        @test @inferred Φ(d,1.0f0) ≈ -3.7462470491110184f-5
        @test Φ(d,1.0f0) isa Float64
        @test @inferred ∇Φ(d,1.0) ≈ 3.38305283586301e-22
        @test ∇Φ(d,1.0) isa Float64
        @test @inferred ∇Φ(d,1.0f0) ≈ 3.38305283586301f-22
        @test ∇Φ(d,1.0f0) isa Float64
        @test @inferred ∇∇Φ(d,1.0) ≈ -2.387252164706999e-22
        @test ∇∇Φ(d,1.0) isa Float64
        @test @inferred ∇∇Φ(d,1.0f0) ≈ -2.387252164706999f-22
        @test ∇∇Φ(d,1.0f0) isa Float64
    end
    @testset "Float32" begin
        @test NFW(1.0f0,1.0f0) isa NFW{Float32}
        @test NFW(1.0f0,1) isa NFW{Float32}
        @test NFW(1,1.0f0) isa NFW{Float32}
        d = NFW(1.0f0,1.0f0)
        @test @inferred ρ(d,1.0) ≈ 0.25f0
        @test ρ(d,1.0) isa Float64
        @test @inferred ρ(d,1.0f0) ≈ 0.25f0
        @test ρ(d,1.0f0) isa Float32
        @test @inferred invρ(d,0.25) ≈ 1.0f0
        @test invρ(d,0.25) isa Float64
        @test @inferred invρ(d,0.25f0) ≈ 1.0f0
        @test invρ(d,0.25f0) isa Float32       
        @test @inferred ∇ρ(d,1.0) ≈ -0.375f0
        @test ∇ρ(d,1.0) isa Float64
        @test @inferred ∇ρ(d,1.0f0) ≈ -0.375f0
        @test ∇ρ(d,1.0f0) isa Float32
        @test @inferred ρmean(d,1.0) ≈ 0.5794415416798359f0
        @test ρmean(d,1.0) isa Float64
        @test @inferred ρmean(d,1.0f0) ≈ 0.5794415416798359f0
        @test ρmean(d,1.0f0) isa Float32
        @test @inferred invρmean(d,0.2) ≈ 1.795108201392415f0
        @test invρmean(d,0.2) isa Float64
        @test @inferred invρmean(d,0.2f0) ≈ 1.795108201392415f0
        @test invρmean(d,0.2f0) isa Float32
        @test @inferred Σ(d,1.0) ≈ Float32(2//3)
        @test Σ(d,1.0) isa Float64
        @test @inferred Σ(d,1.0f0) ≈ Float32(2//3)
        @test Σ(d,1.0f0) isa Float32
        # ∇Σ
        @test @inferred Σmean(d,1.0) ≈ 1.2274112777602189f0
        @test Σmean(d,1.0) isa Float64
        @test @inferred Σmean(d,1.0f0) ≈ 1.2274112777602189f0
        @test Σmean(d,1.0f0) isa Float32
        # ∇Σmean
        @test @inferred invΣ(d,0.5) ≈ 1.259435551338637f0
        @test invΣ(d,0.5) isa Float64
        @test @inferred invΣ(d,0.5f0) ≈ 1.259435551338637f0
        @test invΣ(d,0.5f0) isa Float32
        @test @inferred M(d,1.0) ≈ 2.4271590540348216f0
        @test M(d,1.0) isa Float64
        @test @inferred M(d,1.0f0) ≈ 2.4271590540348216f0
        @test M(d,1.0f0) isa Float32
        @test @inferred ∇M(d,1.0) ≈ Float32(π)
        @test ∇M(d,1.0) isa Float64
        @test @inferred ∇M(d,1.0f0) ≈ Float32(π)
        @test ∇M(d,1.0f0) isa Float32
        @test @inferred invM(d,2.4271590540348216) ≈ 1.0f0
        @test invM(d,2.4271590540348216) isa Float64
        @test @inferred invM(d,2.4271590540348216f0) ≈ 1.0f0
        @test invM(d,2.4271590540348216f0) isa Float32
        # Mtot
        @test @inferred Mproj(d,1.0) ≈ 3.8560262531447647f0
        @test Mproj(d,1.0) isa Float64
        @test @inferred Mproj(d,1.0f0) ≈ 3.8560262531447647f0
        @test Mproj(d,1.0f0) isa Float32
        @test @inferred ∇Mproj(d,1.0) ≈ 4.1887902047863905f0
        @test ∇Mproj(d,1.0) isa Float64
        @test @inferred ∇Mproj(d,1.0f0) ≈ 4.1887902047863905f0
        @test ∇Mproj(d,1.0f0) isa Float32
        @test @inferred invMproj(d,0.2) ≈ 0.11585086460387914f0
        @test invMproj(d,0.2) isa Float64
        @test @inferred invMproj(d,0.2f0) ≈ 0.11585086460387914f0
        @test invMproj(d,0.2f0) isa Float32
        @test @inferred dynamical_time(d,1.0) ≈ 4.753754967831782f11
        @test dynamical_time(d,1.0) isa Float64
        @test @inferred dynamical_time(d,1.0f0) ≈ 4.753754967831782f11
        @test dynamical_time(d,1.0f0) isa Float32
        @test @inferred Vcirc(d,1.0) ≈ 0.003230945727279133f0
        @test Vcirc(d,1.0) isa Float64
        @test @inferred Vcirc(d,1.0f0) ≈ 0.003230945727279133f0
        @test Vcirc(d,1.0f0) isa Float32
        @test @inferred Vesc(d,1.0) ≈ 0.008655919418653363f0
        @test Vesc(d,1.0) isa Float64
        @test @inferred Vesc(d,1.0f0) ≈ 0.008655919418653363f0
        @test Vesc(d,1.0f0) isa Float32
        @test @inferred all(Vmax(d) .≈ (0.003418455956363109f0, 2.1625815870646097f0))
        @test Vmax(d) isa NTuple{2,Float32}
        @test @inferred Φ(d,1.0) ≈ -3.7462470491110184f-5
        @test Φ(d,1.0) isa Float64
        @test @inferred Φ(d,1.0f0) ≈ -3.7462470491110184f-5
        @test Φ(d,1.0f0) isa Float32
        @test @inferred ∇Φ(d,1.0) ≈ 3.38305283586301f-22
        @test ∇Φ(d,1.0) isa Float64
        @test @inferred ∇Φ(d,1.0f0) ≈ 3.38305283586301f-22
        @test ∇Φ(d,1.0f0) isa Float32
        @test @inferred ∇∇Φ(d,1.0) ≈ -2.387252164706999f-22
        @test ∇∇Φ(d,1.0) isa Float64
        @test @inferred ∇∇Φ(d,1.0f0) ≈ -2.387252164706999f-22
        @test ∇∇Φ(d,1.0f0) isa Float32
    end
end
