# this testset will be exhaustive; if the units work on this one, they should work on most other ones too. 
@testset "GeneralIsothermal Units" begin
    @testset "Float64" begin
        @test GeneralIsothermal(1.0*ua.Msun/ua.kpc^3,1.0*ua.kpc,1.0) isa GeneralIsothermal{Float64}
        @test GeneralIsothermal(1.0*ua.Msun/ua.kpc^3,1.0*ua.kpc,1) isa GeneralIsothermal{Float64}
        @test GeneralIsothermal(1.0*ua.Msun/ua.kpc^3,1*ua.kpc,1.0f0) isa GeneralIsothermal{Float64}
        @test GeneralIsothermal(1.0*u.kg/u.m^3,1.0*u.m,1.0) isa GeneralIsothermal{Float64}
        # this should construct d = GeneralIsothermal(1.0,1.0,1.0); then run tests with Reals
        d = GeneralIsothermal(1.0*GalaxyProfiles.defaultunits.density,1.0*GalaxyProfiles.defaultunits.length,1.0)
        @test d == GeneralIsothermal(1.0,1.0,1.0)
        @test @inferred ρ(d,1.0) == 1.0 
        @test ρ(d,1.0) isa Float64
        @test @inferred ρ(d,1.0f0) == 1.0
        @test ρ(d,1.0f0) isa Float64
        @test @inferred invρ(d,1.0) == 1.0
        @test invρ(d,1.0) isa Float64
        @test @inferred invρ(d,1.0f0) == 1.0
        @test invρ(d,1.0f0) isa Float64
        @test @inferred ∇ρ(d,1.0) == -1.0
        @test ∇ρ(d,1.0) isa Float64
        @test @inferred ∇ρ(d,1.0f0) == -1.0
        @test ∇ρ(d,1.0f0) isa Float64
        @test @inferred Σ(d,1.0) == Inf
        @test Σ(d,1.0) isa Float64
        @test @inferred Σ(d,1.0f0) == Inf
        @test Σ(d,1.0f0) isa Float64
        # now run tests with Unitful inputs
        @test @inferred ρ(GalaxyProfiles.defaultunits.density,d,1.0) == 1.0*GalaxyProfiles.defaultunits.density
        @test ρ(GalaxyProfiles.defaultunits.density,d,1.0) |> u.ustrip isa Float64
        @test ρ(GalaxyProfiles.defaultunits.density,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred ρ(GalaxyProfiles.defaultunits.density,d,1.0*GalaxyProfiles.defaultunits.length) == 1.0*GalaxyProfiles.defaultunits.density
        @test ρ(GalaxyProfiles.defaultunits.density,d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test ρ(GalaxyProfiles.defaultunits.density,d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test @inferred ρ(d,1.0*GalaxyProfiles.defaultunits.length) == 1.0*GalaxyProfiles.defaultunits.density
        @test ρ(d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test ρ(d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        ####################################
        @test @inferred invρ(GalaxyProfiles.defaultunits.length,d,1.0) == 1.0*GalaxyProfiles.defaultunits.length
        @test invρ(GalaxyProfiles.defaultunits.length,d,1.0) |> u.ustrip isa Float64
        @test invρ(GalaxyProfiles.defaultunits.length,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred invρ(GalaxyProfiles.defaultunits.length,d,1.0*GalaxyProfiles.defaultunits.density) == 1.0*GalaxyProfiles.defaultunits.length
        @test invρ(GalaxyProfiles.defaultunits.length,d,1.0*GalaxyProfiles.defaultunits.density) |> u.ustrip isa Float64
        @test invρ(GalaxyProfiles.defaultunits.length,d,1.0f0*GalaxyProfiles.defaultunits.density) |> u.ustrip isa Float64
        @test @inferred invρ(d,1.0*GalaxyProfiles.defaultunits.density) == 1.0*GalaxyProfiles.defaultunits.length
        @test invρ(d,1.0*GalaxyProfiles.defaultunits.density) |> u.ustrip isa Float64
        @test invρ(d,1.0f0*GalaxyProfiles.defaultunits.density) |> u.ustrip isa Float64
        ####################################
        @test @inferred ∇ρ(GalaxyProfiles.defaultunits.∇density,d,1.0) == -1.0*GalaxyProfiles.defaultunits.∇density
        @test ∇ρ(GalaxyProfiles.defaultunits.∇density,d,1.0) |> u.ustrip isa Float64
        @test ∇ρ(GalaxyProfiles.defaultunits.∇density,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred ∇ρ(GalaxyProfiles.defaultunits.∇density,d,1.0*GalaxyProfiles.defaultunits.length) == -1.0*GalaxyProfiles.defaultunits.∇density
        @test ∇ρ(GalaxyProfiles.defaultunits.∇density,d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test ∇ρ(GalaxyProfiles.defaultunits.∇density,d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test @inferred ∇ρ(d,1.0*GalaxyProfiles.defaultunits.length) == -1.0*GalaxyProfiles.defaultunits.∇density
        @test ∇ρ(d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test ∇ρ(d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        ####################################
        @test @inferred Σ(GalaxyProfiles.defaultunits.surfacedensity,d,1.0) == Inf*GalaxyProfiles.defaultunits.surfacedensity
        @test Σ(GalaxyProfiles.defaultunits.surfacedensity,d,1.0) |> u.ustrip isa Float64
        @test Σ(GalaxyProfiles.defaultunits.surfacedensity,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred Σ(GalaxyProfiles.defaultunits.surfacedensity,d,1.0*GalaxyProfiles.defaultunits.length) == Inf*GalaxyProfiles.defaultunits.surfacedensity
        @test Σ(GalaxyProfiles.defaultunits.surfacedensity,d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test Σ(GalaxyProfiles.defaultunits.surfacedensity,d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test @inferred Σ(d,1.0*GalaxyProfiles.defaultunits.length) == Inf*GalaxyProfiles.defaultunits.surfacedensity
        @test Σ(d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64        
        @test Σ(d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64        

        # this should construct d = GeneralIsothermal(1.0,1.0,2.0); then run tests with Reals
        d = GeneralIsothermal(1.0*GalaxyProfiles.defaultunits.density,1.0*GalaxyProfiles.defaultunits.length,2.0)
        @test d == GeneralIsothermal(1.0,1.0,2.0)
        
        @test @inferred invΣ(d,π) == 1.0
        @test invΣ(d,π) isa Float64
        @test @inferred invΣ(d,Float32(π)) ≈ 1.0f0
        @test invΣ(d,Float32(π)) isa Float64
        @test @inferred ∇Σ(d,1.0) == -Float64(π)
        @test ∇Σ(d,1.0) isa Float64
        @test @inferred ∇Σ(d,1.0f0) == -Float64(π)
        @test ∇Σ(d,1.0f0) isa Float64
        @test @inferred M(d,1.0) == 4π
        @test M(d,1.0) isa Float64
        @test @inferred M(d,1.0f0) == 4π
        @test M(d,1.0f0) isa Float64
        @test @inferred invM(d,4π) == 1.0
        @test invM(d,4π) isa Float64
        @test @inferred invM(d,Float32(4π)) ≈ 1.0f0
        @test invM(d,Float32(4π)) isa Float64
        @test @inferred ∇M(d,1.0) == 4π
        @test ∇M(d,1.0) isa Float64
        @test @inferred ∇M(d,1.0f0) == 4π
        @test ∇M(d,1.0f0) isa Float64
        @test @inferred Mproj(d,1.0) == 2*π^2
        @test Mproj(d,1.0) isa Float64
        @test @inferred Mproj(d,1.0f0) == 2*π^2
        @test Mproj(d,1.0f0) isa Float64
        @test @inferred ∇Mproj(d,1.0) == 2*π^2
        @test ∇Mproj(d,1.0) isa Float64
        @test @inferred ∇Mproj(d,1.0f0) == 2*π^2
        @test ∇Mproj(d,1.0f0) isa Float64
        @test @inferred invMproj(d,2*π^2) == 1.0
        @test invMproj(d,2*π^2) isa Float64
        @test @inferred invMproj(d,Float32(2*π^2)) ≈ 1.0f0
        @test invMproj(d,Float32(2*π^2)) isa Float64
        # @test @inferred Vcirc(d,1.0) == 1.0
        @test Vcirc(d,1.0) isa Float64
        # @test @inferred Vcirc(d,1.0) ≈ 1.0f0
        @test Vcirc(d,Float32(1.0)) isa Float64
        @test_throws DomainError Vesc(d,1.0)   # α=2 is not valid for Vesc
        let d = GeneralIsothermal(1.0,1.0,2.5) # introduce new d just for Vesc
            # @test @inferred Vesc(d,1.0) == 1.0
            @test Vesc(d,1.0) isa Float64
            # @test @inferred Vesc(d,1.0) ≈ 1.0f0
            @test Vesc(d,Float32(1.0)) isa Float64
        end
        # @test @inferred Φ(d,1.0) == 1.0
        @test Φ(d,1.0) isa Float64
        # @test @inferred Φ(d,1.0) ≈ 1.0f0
        @test Φ(d,Float32(1.0)) isa Float64
        # @test @inferred ∇Φ(d,1.0) == 1.0
        @test ∇Φ(d,1.0) isa Float64
        # @test @inferred ∇Φ(d,1.0) ≈ 1.0f0
        @test ∇Φ(d,Float32(1.0)) isa Float64
        # @test @inferred ∇∇Φ(d,1.0) == 1.0
        @test ∇∇Φ(d,1.0) isa Float64
        # @test @inferred ∇∇Φ(d,1.0) ≈ 1.0f0
        @test ∇∇Φ(d,Float32(1.0)) isa Float64

        # Now run the above tests with Unitful inputs
        @test @inferred invΣ(GalaxyProfiles.defaultunits.length,d,π) == 1.0*GalaxyProfiles.defaultunits.length
        @test invΣ(GalaxyProfiles.defaultunits.length,d,π) |> u.ustrip isa Float64
        @test invΣ(GalaxyProfiles.defaultunits.length,d,Float32(π)) |> u.ustrip isa Float64
        @test @inferred invΣ(GalaxyProfiles.defaultunits.length,d,π*GalaxyProfiles.defaultunits.surfacedensity) == 1*GalaxyProfiles.defaultunits.length
        @test invΣ(GalaxyProfiles.defaultunits.length,d,π*GalaxyProfiles.defaultunits.surfacedensity) |> u.ustrip isa Float64
        @test invΣ(GalaxyProfiles.defaultunits.length,d,Float32(π)*GalaxyProfiles.defaultunits.surfacedensity) |> u.ustrip isa Float64
        @test @inferred invΣ(d,π*GalaxyProfiles.defaultunits.surfacedensity) == 1.0*GalaxyProfiles.defaultunits.length
        @test invΣ(d,π*GalaxyProfiles.defaultunits.surfacedensity) |> u.ustrip isa Float64 
        @test invΣ(d,Float32(π)*GalaxyProfiles.defaultunits.surfacedensity) |> u.ustrip isa Float64
        ####################################
        @test @inferred ∇Σ(GalaxyProfiles.defaultunits.density,d,1.0) == -π*GalaxyProfiles.defaultunits.density
        @test ∇Σ(GalaxyProfiles.defaultunits.density,d,1.0) |> u.ustrip isa Float64
        @test ∇Σ(GalaxyProfiles.defaultunits.density,d,Float32(1.0)) |> u.ustrip isa Float64
        @test @inferred ∇Σ(GalaxyProfiles.defaultunits.density,d,1.0*GalaxyProfiles.defaultunits.length) == -π*GalaxyProfiles.defaultunits.density
        @test ∇Σ(GalaxyProfiles.defaultunits.density,d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test ∇Σ(GalaxyProfiles.defaultunits.density,d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test @inferred ∇Σ(d,1.0*GalaxyProfiles.defaultunits.length) == -π*GalaxyProfiles.defaultunits.density
        @test ∇Σ(d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64 
        @test ∇Σ(d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        ####################################
        @test @inferred M(GalaxyProfiles.defaultunits.mass,d,1.0) == 4π*GalaxyProfiles.defaultunits.mass
        @test M(GalaxyProfiles.defaultunits.mass,d,1.0) |> u.ustrip isa Float64
        @test M(GalaxyProfiles.defaultunits.mass,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred M(GalaxyProfiles.defaultunits.mass,d,1.0*GalaxyProfiles.defaultunits.length) == 4π*GalaxyProfiles.defaultunits.mass
        @test M(GalaxyProfiles.defaultunits.mass,d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test M(GalaxyProfiles.defaultunits.mass,d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test @inferred M(d,1.0*GalaxyProfiles.defaultunits.length) == 4π*GalaxyProfiles.defaultunits.mass
        @test M(d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test M(d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        ####################################
        @test @inferred invM(GalaxyProfiles.defaultunits.length,d,4π) == 1.0*GalaxyProfiles.defaultunits.length
        @test invM(GalaxyProfiles.defaultunits.length,d,4π) |> u.ustrip isa Float64
        @test invM(GalaxyProfiles.defaultunits.length,d,Float32(4π)) |> u.ustrip isa Float64
        @test @inferred invM(GalaxyProfiles.defaultunits.length,d,4π*GalaxyProfiles.defaultunits.mass) == 1.0*GalaxyProfiles.defaultunits.length
        @test invM(GalaxyProfiles.defaultunits.length,d,4π*GalaxyProfiles.defaultunits.mass) |> u.ustrip isa Float64
        @test invM(GalaxyProfiles.defaultunits.length,d,Float32(4π)*GalaxyProfiles.defaultunits.mass) |> u.ustrip isa Float64
        @test @inferred invM(d,4π*GalaxyProfiles.defaultunits.mass) == 1.0*GalaxyProfiles.defaultunits.length
        @test invM(d,4π*GalaxyProfiles.defaultunits.mass) |> u.ustrip isa Float64
        @test invM(d,Float32(4π)*GalaxyProfiles.defaultunits.mass) |> u.ustrip isa Float64
        ####################################
        @test @inferred ∇M(GalaxyProfiles.defaultunits.∇mass,d,1.0) == 4π*GalaxyProfiles.defaultunits.∇mass
        @test ∇M(GalaxyProfiles.defaultunits.∇mass,d,1.0) |> u.ustrip isa Float64
        @test ∇M(GalaxyProfiles.defaultunits.∇mass,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred ∇M(GalaxyProfiles.defaultunits.∇mass,d,1.0*GalaxyProfiles.defaultunits.length) == 4π*GalaxyProfiles.defaultunits.∇mass
        @test ∇M(GalaxyProfiles.defaultunits.∇mass,d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test ∇M(GalaxyProfiles.defaultunits.∇mass,d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test @inferred ∇M(d,1.0*GalaxyProfiles.defaultunits.length) == 4π*GalaxyProfiles.defaultunits.∇mass
        @test ∇M(d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test ∇M(d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        ####################################
        @test @inferred Mproj(GalaxyProfiles.defaultunits.mass,d,1.0) == 2*π^2 * GalaxyProfiles.defaultunits.mass
        @test Mproj(GalaxyProfiles.defaultunits.mass,d,1.0) |> u.ustrip isa Float64
        @test Mproj(GalaxyProfiles.defaultunits.mass,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred Mproj(GalaxyProfiles.defaultunits.mass,d,1.0*GalaxyProfiles.defaultunits.length) == 2*π^2 * GalaxyProfiles.defaultunits.mass
        @test Mproj(GalaxyProfiles.defaultunits.mass,d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test Mproj(GalaxyProfiles.defaultunits.mass,d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test @inferred Mproj(d,1.0*GalaxyProfiles.defaultunits.length) == 2*π^2 * GalaxyProfiles.defaultunits.mass
        @test Mproj(d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test Mproj(d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        ####################################
        @test @inferred ∇Mproj(GalaxyProfiles.defaultunits.∇mass,d,1.0) == 2*π^2*GalaxyProfiles.defaultunits.∇mass
        @test ∇Mproj(GalaxyProfiles.defaultunits.∇mass,d,1.0) |> u.ustrip isa Float64
        @test ∇Mproj(GalaxyProfiles.defaultunits.∇mass,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred ∇Mproj(GalaxyProfiles.defaultunits.∇mass,d,1.0*GalaxyProfiles.defaultunits.length) == 2*π^2*GalaxyProfiles.defaultunits.∇mass
        @test ∇Mproj(GalaxyProfiles.defaultunits.∇mass,d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test ∇Mproj(GalaxyProfiles.defaultunits.∇mass,d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test @inferred ∇Mproj(d,1.0*GalaxyProfiles.defaultunits.length) == 2*π^2*GalaxyProfiles.defaultunits.∇mass
        @test ∇Mproj(d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test ∇Mproj(d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        ####################################
        @test @inferred invMproj(GalaxyProfiles.defaultunits.length,d,2*π^2) == 1.0*GalaxyProfiles.defaultunits.length
        @test invMproj(GalaxyProfiles.defaultunits.length,d,2*π^2) |> u.ustrip isa Float64
        @test invMproj(GalaxyProfiles.defaultunits.length,d,Float32(2*π^2)) |> u.ustrip isa Float64
        @test @inferred invMproj(GalaxyProfiles.defaultunits.length,d,2*π^2*GalaxyProfiles.defaultunits.mass) == 1.0*GalaxyProfiles.defaultunits.length
        @test invMproj(GalaxyProfiles.defaultunits.length,d,2*π^2*GalaxyProfiles.defaultunits.mass) |> u.ustrip isa Float64
        @test invMproj(GalaxyProfiles.defaultunits.length,d,Float32(2*π^2)*GalaxyProfiles.defaultunits.mass) |> u.ustrip isa Float64
        @test @inferred invMproj(d,2*π^2*GalaxyProfiles.defaultunits.mass) == 1.0*GalaxyProfiles.defaultunits.length
        @test invMproj(d,2*π^2*GalaxyProfiles.defaultunits.mass) |> u.ustrip isa Float64
        @test invMproj(d,Float32(2*π^2)*GalaxyProfiles.defaultunits.mass) |> u.ustrip isa Float64
        ####################################
        @test @inferred Vcirc(GalaxyProfiles.defaultunits.velocity,d,1.0) == Vcirc(d,1.0) * GalaxyProfiles.defaultunits.velocity
        @test Vcirc(GalaxyProfiles.defaultunits.velocity,d,1.0) |> u.ustrip isa Float64
        @test Vcirc(GalaxyProfiles.defaultunits.velocity,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred Vcirc(GalaxyProfiles.defaultunits.velocity,d,1.0*GalaxyProfiles.defaultunits.length) == Vcirc(d,1.0) * GalaxyProfiles.defaultunits.velocity
        @test Vcirc(GalaxyProfiles.defaultunits.velocity,d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test Vcirc(GalaxyProfiles.defaultunits.velocity,d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test @inferred Vcirc(d,1.0*GalaxyProfiles.defaultunits.length) == Vcirc(d,1.0) * GalaxyProfiles.defaultunits.velocity
        @test Vcirc(d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test Vcirc(d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        ####################################
        @test_throws DomainError Vesc(GalaxyProfiles.defaultunits.velocity,d,1.0)   # α=2 is not valid for Vesc
        @test_throws DomainError Vesc(GalaxyProfiles.defaultunits.velocity,d,1.0*GalaxyProfiles.defaultunits.length)   # α=2 is not valid for Vesc
        @test_throws DomainError Vesc(d,1.0*GalaxyProfiles.defaultunits.length)   # α=2 is not valid for Vesc
        ####################################
        let d = GeneralIsothermal(1.0*GalaxyProfiles.defaultunits.density,1.0*GalaxyProfiles.defaultunits.length,2.5) # introduce new d just for Vesc
            @test d == GeneralIsothermal(1.0,1.0,2.5)
            @test @inferred Vesc(GalaxyProfiles.defaultunits.velocity,d,1.0) == Vesc(d,1.0) * GalaxyProfiles.defaultunits.velocity
            @test Vesc(GalaxyProfiles.defaultunits.velocity,d,1.0) |> u.ustrip isa Float64
            @test Vesc(GalaxyProfiles.defaultunits.velocity,d,1.0f0) |> u.ustrip isa Float64
            @test @inferred Vesc(GalaxyProfiles.defaultunits.velocity,d,1.0*GalaxyProfiles.defaultunits.length) == Vesc(d,1.0) * GalaxyProfiles.defaultunits.velocity
            @test Vesc(GalaxyProfiles.defaultunits.velocity,d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
            @test Vesc(GalaxyProfiles.defaultunits.velocity,d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
            @test @inferred Vesc(d,1.0*GalaxyProfiles.defaultunits.length) == Vesc(d,1.0) * GalaxyProfiles.defaultunits.velocity
            @test Vesc(d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
            @test Vesc(d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        end
        ####################################
        # begin # potentials; for some reason these mess up Emacs julia mode so stick them in a begin block
        @test @inferred Φ(GalaxyProfiles.defaultunits.Φunit,d,1.0) == Φ(d,1.0) * GalaxyProfiles.defaultunits.Φunit
        @test Φ(GalaxyProfiles.defaultunits.Φunit,d,1.0) |> u.ustrip isa Float64
        @test Φ(GalaxyProfiles.defaultunits.Φunit,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred Φ(GalaxyProfiles.defaultunits.Φunit,d,1.0*GalaxyProfiles.defaultunits.length) == Φ(d,1.0) * GalaxyProfiles.defaultunits.Φunit
        @test Φ(GalaxyProfiles.defaultunits.Φunit,d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test Φ(GalaxyProfiles.defaultunits.Φunit,d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test @inferred Φ(d,1.0*GalaxyProfiles.defaultunits.length) == Φ(d,1.0) * GalaxyProfiles.defaultunits.Φunit
        @test Φ(d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test Φ(d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        ####################################
        @test @inferred ∇Φ(GalaxyProfiles.defaultunits.∇Φunit,d,1.0) == ∇Φ(d,1.0) * GalaxyProfiles.defaultunits.∇Φunit
        @test @inferred ∇Φ(GalaxyProfiles.defaultunits.∇Φunit,d,1.0) == ∇Φ(d,1.0) * GalaxyProfiles.defaultunits.∇Φunit
        @test ∇Φ(GalaxyProfiles.defaultunits.∇Φunit,d,1.0) |> u.ustrip isa Float64
        @test ∇Φ(GalaxyProfiles.defaultunits.∇Φunit,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred ∇Φ(GalaxyProfiles.defaultunits.∇Φunit,d,1.0*GalaxyProfiles.defaultunits.length) == ∇Φ(d,1.0) * GalaxyProfiles.defaultunits.∇Φunit
        @test ∇Φ(GalaxyProfiles.defaultunits.∇Φunit,d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test ∇Φ(GalaxyProfiles.defaultunits.∇Φunit,d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test @inferred ∇Φ(d,1.0*GalaxyProfiles.defaultunits.length) == ∇Φ(d,1.0) * GalaxyProfiles.defaultunits.∇Φunit
        @test ∇Φ(d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test ∇Φ(d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        ####################################
        @test @inferred ∇∇Φ(GalaxyProfiles.defaultunits.∇∇Φunit,d,1.0) == ∇∇Φ(d,1.0) * GalaxyProfiles.defaultunits.∇∇Φunit
        @test ∇∇Φ(GalaxyProfiles.defaultunits.∇∇Φunit,d,1.0) |> u.ustrip isa Float64
        @test ∇∇Φ(GalaxyProfiles.defaultunits.∇∇Φunit,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred ∇∇Φ(GalaxyProfiles.defaultunits.∇∇Φunit,d,1.0*GalaxyProfiles.defaultunits.length) == ∇∇Φ(d,1.0) * GalaxyProfiles.defaultunits.∇∇Φunit
        @test ∇∇Φ(GalaxyProfiles.defaultunits.∇∇Φunit,d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test ∇∇Φ(GalaxyProfiles.defaultunits.∇∇Φunit,d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test @inferred ∇∇Φ(d,1.0*GalaxyProfiles.defaultunits.length) == ∇∇Φ(d,1.0) * GalaxyProfiles.defaultunits.∇∇Φunit
        @test ∇∇Φ(d,1.0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        @test ∇∇Φ(d,1.0f0*GalaxyProfiles.defaultunits.length) |> u.ustrip isa Float64
        ####################################
    end

    # float32 behavior has already been validated for the Real components; just do a few tests here
    @testset "Float32" begin
        @test GeneralIsothermal(1.0f0 * GalaxyProfiles.defaultunits.density,1.0f0 * GalaxyProfiles.defaultunits.length,1.0f0) isa GeneralIsothermal{Float32}
        @test GeneralIsothermal(1.0f0 * GalaxyProfiles.defaultunits.density,1.0f0 * GalaxyProfiles.defaultunits.length,1) isa GeneralIsothermal{Float32}
    end
end
