using GalaxyProfiles
using Test

@testset "ExponentialDisk" begin
    @testset "Float64" begin
        @test ExponentialDisk(1.0,1.0) isa ExponentialDisk{Float64}
        @test ExponentialDisk(1.0,1) isa ExponentialDisk{Float64}
        @test ExponentialDisk(1.0,1.0f0) isa ExponentialDisk{Float64}
        @test ExponentialDisk(1.0; M=1.0) isa ExponentialDisk{Float64}
        @test ExponentialDisk(1.0; Σ0=1.0) isa ExponentialDisk{Float64}
        @test ExponentialDiskDHI(1.0,1.0) isa ExponentialDisk{Float64}
        @test ExponentialDiskDHI(1.0,1.0,1e6) isa ExponentialDisk{Float64}
        @test ExponentialDiskDHI(1.0,1.0f0,1e6) isa ExponentialDisk{Float64}
        d = ExponentialDisk(1.0,1.0)
        @test @inferred Σ(d,1.0) == exp(-1)
        @test Σ(d,1.0) isa Float64
        @test @inferred Σ(d,1.0f0) == exp(-1)
        @test Σ(d,1.0f0) isa Float64
        @test @inferred invΣ(d,exp(-1)) == 1.0
        @test invΣ(d,exp(-1)) isa Float64
        @test @inferred invΣ(d,Float32(exp(-1))) ≈ 1.0f0
        @test invΣ(d,Float32(exp(-1))) isa Float64
        @test @inferred ∇Σ(d,1.0) == -exp(-1)
        @test ∇Σ(d,1.0) isa Float64
        @test @inferred ∇Σ(d,1.0f0) == -exp(-1)
        @test ∇Σ(d,1.0f0) isa Float64
        @test @inferred Mproj(d,1.0) == 2π * (1 - 2*exp(-1))
        @test Mproj(d,1.0) isa Float64
        @test @inferred Mproj(d,1.0f0) == 2π * (1 - 2*exp(-1))
        @test Mproj(d,1.0f0) isa Float64
        @test @inferred ∇Mproj(d,1.0) == 2π * exp(-1)
        @test ∇Mproj(d,1.0) isa Float64
        @test @inferred ∇Mproj(d,1.0f0) == 2π * exp(-1)
        @test ∇Mproj(d,1.0f0) isa Float64
        @test @inferred invMproj(d,2π * (1 - 2*exp(-1))) ≈ 1.0
        @test invMproj(d,2π * (1 - 2*exp(-1))) isa Float64
        @test @inferred invMproj(d,Float32(2π * (1 - 2*exp(-1)))) ≈ 1.0f0
        @test invMproj(d,Float32(2π * (1 - 2*exp(-1)))) isa Float64
        @test @inferred Mtot(d) == 2π
        @test Mtot(d) isa Float64
    end
    @testset "Float32" begin
        @test ExponentialDisk(1.0f0,1.0f0) isa ExponentialDisk{Float32}
        @test ExponentialDisk(1.0f0,1) isa ExponentialDisk{Float32}
        @test ExponentialDisk(1,1.0f0) isa ExponentialDisk{Float32}
        @test ExponentialDisk(1.0f0; M=1.0f0) isa ExponentialDisk{Float32}
        @test ExponentialDisk(1.0f0; Σ0=1.0f0) isa ExponentialDisk{Float32}
        @test ExponentialDiskDHI(1.0f0,1.0f0) isa ExponentialDisk{Float32}
        @test ExponentialDiskDHI(1.0f0,1.0f0,1f6) isa ExponentialDisk{Float32}
        d = ExponentialDisk(1.0f0,1.0f0)
        @test @inferred Σ(d,1.0) == exp(-1)
        @test Σ(d,1.0) isa Float64
        @test @inferred Σ(d,1.0f0) == Float32(exp(-1))
        @test Σ(d,1.0f0) isa Float32
        @test @inferred invΣ(d,exp(-1)) == 1.0
        @test invΣ(d,exp(-1)) isa Float64
        @test @inferred invΣ(d,Float32(exp(-1))) ≈ 1.0f0
        @test invΣ(d,Float32(exp(-1))) isa Float32
        @test @inferred ∇Σ(d,1.0) == -exp(-1)
        @test ∇Σ(d,1.0) isa Float64
        @test @inferred ∇Σ(d,1.0f0) == Float32(-exp(-1))
        @test ∇Σ(d,1.0f0) isa Float32
        @test @inferred Mproj(d,1.0) == 2π * (1 - 2*exp(-1))
        @test Mproj(d,1.0) isa Float64
        @test @inferred Mproj(d,1.0f0) ≈ Float32(2π * (1 - 2*exp(-1)))
        @test Mproj(d,1.0f0) isa Float32
        @test @inferred ∇Mproj(d,1.0) == 2π * exp(-1)
        @test ∇Mproj(d,1.0) isa Float64
        @test @inferred ∇Mproj(d,1.0f0) == Float32(2π * exp(-1))
        @test ∇Mproj(d,1.0f0) isa Float32
        @test @inferred invMproj(d,2π * (1 - 2*exp(-1))) ≈ 1.0f0
        @test invMproj(d,2π * (1 - 2*exp(-1))) isa Float64
        @test @inferred invMproj(d,Float32(2π * (1 - 2*exp(-1)))) ≈ 1.0f0
        @test invMproj(d,Float32(2π * (1 - 2*exp(-1)))) isa Float32
        @test @inferred Mtot(d) == Float32(2π)
        @test Mtot(d) isa Float32
    end
end
