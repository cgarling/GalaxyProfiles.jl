@testset "GeneralIsothermal" begin
    @testset "Float64" begin
        @test GeneralIsothermal(1.0,1.0,1.0) isa GeneralIsothermal{Float64}
        @test GeneralIsothermal(1.0,1.0,1) isa GeneralIsothermal{Float64}
        @test GeneralIsothermal(1.0,1.0,1.0f0) isa GeneralIsothermal{Float64}
        d = GeneralIsothermal(1.0,1.0,1.0)
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
        
        d = GeneralIsothermal(1.0,1.0,2.0)
        
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
    end
    @testset "Float32" begin
        @test GeneralIsothermal(1.0f0,1.0f0,1.0f0) isa GeneralIsothermal{Float32}
        @test GeneralIsothermal(1.0f0,1.0f0,1) isa GeneralIsothermal{Float32}
        d = GeneralIsothermal(1.0f0,1.0f0,1.0f0)
        @test @inferred ρ(d,1.0) == 1.0
        @test ρ(d,1.0) isa Float64
        @test @inferred ρ(d,1.0f0) == 1.0f0
        @test ρ(d,1.0f0) isa Float32
        @test @inferred invρ(d,1.0) == 1.0
        @test invρ(d,1.0) isa Float64
        @test @inferred invρ(d,1.0f0) == 1.0f0
        @test invρ(d,1.0f0) isa Float32
        @test @inferred ∇ρ(d,1.0) == -1.0
        @test ∇ρ(d,1.0) isa Float64
        @test @inferred ∇ρ(d,1.0f0) == -1.0f0
        @test ∇ρ(d,1.0f0) isa Float32
        @test @inferred Σ(d,1.0) == Inf
        @test Σ(d,1.0) isa Float64
        @test @inferred Σ(d,1.0f0) == Inf32
        @test Σ(d,1.0f0) isa Float32
        
        d = GeneralIsothermal(1.0f0,1.0f0,2.0f0)
        
        @test @inferred invΣ(d,π) == 1.0f0
        @test invΣ(d,π) isa Float32
        @test @inferred invΣ(d,Float64(π)) == 1.0
        @test invΣ(d,Float64(π)) isa Float64
        @test @inferred invΣ(d,Float32(π)) == 1.0f0
        @test invΣ(d,Float32(π)) isa Float32
        @test @inferred ∇Σ(d,1.0) == -Float64(π)
        @test ∇Σ(d,1.0) isa Float64
        @test @inferred ∇Σ(d,1.0f0) == -Float32(π)
        @test ∇Σ(d,1.0f0) isa Float32
        @test @inferred M(d,1.0) == 4π
        @test M(d,1.0) isa Float64
        @test @inferred M(d,1.0f0) == Float32(4π)
        @test M(d,1.0f0) isa Float32
        @test @inferred invM(d,4π) == 1.0
        @test invM(d,4π) isa Float64
        @test @inferred invM(d,Float32(4π)) == 1.0f0
        @test invM(d,Float32(4π)) isa Float32
        @test @inferred ∇M(d,1.0) == 4π
        @test ∇M(d,1.0) isa Float64
        @test @inferred ∇M(d,1.0f0) == Float32(4π)
        @test ∇M(d,1.0f0) isa Float32
        @test @inferred Mproj(d,1.0) == 2*π^2
        @test Mproj(d,1.0) isa Float64
        @test @inferred Mproj(d,1.0f0) == Float32(2*π^2)
        @test Mproj(d,1.0f0) isa Float32
        @test @inferred ∇Mproj(d,1.0) == 2*π^2
        @test ∇Mproj(d,1.0) isa Float64
        @test @inferred ∇Mproj(d,1.0f0) == Float32(2*π^2)
        @test ∇Mproj(d,1.0f0) isa Float32
        @test @inferred invMproj(d,2*π^2) == 1.0
        @test invMproj(d,2*π^2) isa Float64
        @test @inferred invMproj(d,Float32(2*π^2)) == 1.0f0
        @test invMproj(d,Float32(2*π^2)) isa Float32
        # @test @inferred Vcirc(d,1.0) == 1.0
        @test Vcirc(d,1.0) isa Float64
        # @test @inferred Vcirc(d,1.0) ≈ 1.0f0
        @test Vcirc(d,Float32(1.0)) isa Float32
        @test_throws DomainError Vesc(d,1.0)   # α=2 is not valid for Vesc
        let d = GeneralIsothermal(1.0f0,1.0f0,2.5f0) # introduce new d just for Vesc
            # @test @inferred Vesc(d,1.0) == 1.0
            @test Vesc(d,1.0) isa Float64
            # @test @inferred Vesc(d,1.0) ≈ 1.0f0
            @test Vesc(d,Float32(1.0)) isa Float32
        end
        # @test @inferred Φ(d,1.0) == 1.0
        @test Φ(d,1.0) isa Float64
        # @test @inferred Φ(d,1.0) ≈ 1.0f0
        @test Φ(d,Float32(1.0)) isa Float32
        # @test @inferred ∇Φ(d,1.0) == 1.0
        @test ∇Φ(d,1.0) isa Float64
        # @test @inferred ∇Φ(d,1.0) ≈ 1.0f0
        @test ∇Φ(d,Float32(1.0)) isa Float32
        # @test @inferred ∇∇Φ(d,1.0) == 1.0
        @test ∇∇Φ(d,1.0) isa Float64
        # @test @inferred ∇∇Φ(d,1.0) ≈ 1.0f0
        @test ∇∇Φ(d,Float32(1.0)) isa Float32
    end
end
