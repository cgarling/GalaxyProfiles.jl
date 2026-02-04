# Test type constructors and common methods with Unitful.jl arguments

using GalaxyProfiles
using Test
import Unitful as u
import UnitfulAstro as ua

if isdefined(Base, :get_extension)
    ext = Base.get_extension(GalaxyProfiles, :GalaxyProfilesUnitfulExt)
    defaultunits = ext.defaultunits
else # For Julia < 1.9 without package extensions
    defaultunits = GalaxyProfiles.GalaxyProfilesUnitfulExt.defaultunits
end

@testset "Unitful Constructors" begin
    @testset "NFW" begin
        @testset "Float64" begin
            # Units are now PRESERVED in struct fields
            nfw_unitful = NFW(1.0*defaultunits.density, 1.0*defaultunits.length)
            @test nfw_unitful.ρ0 isa u.Quantity{Float64}
            @test nfw_unitful.rs isa u.Quantity{Float64}
            # Plain number constructor still works
            nfw_plain = NFW(1.0, 1.0)
            @test nfw_plain.ρ0 == 1.0
            @test nfw_plain.rs == 1.0
        end
        @testset "Float32" begin
            nfw_unitful = NFW(1.0f0*defaultunits.density, 1.0f0*defaultunits.length)
            @test nfw_unitful.ρ0 isa u.Quantity{Float32}
            @test nfw_unitful.rs isa u.Quantity{Float32}
        end
    end
    
    @testset "ExponentialDisk" begin
        @testset "Float64" begin
            # Units are now PRESERVED in struct fields
            ed_unitful = ExponentialDisk(1.0*defaultunits.surfacedensity, 1.0*defaultunits.length)
            @test ed_unitful.Σ0 isa u.Quantity{Float64}
            @test ed_unitful.rs isa u.Quantity{Float64}
            # Plain number constructor still works
            ed_plain = ExponentialDisk(1.0, 1.0)
            @test ed_plain.Σ0 == 1.0
            @test ed_plain.rs == 1.0
            # Test keyword constructors with Unitful
            ed_M = ExponentialDisk(1.0*defaultunits.length; M=1.0*defaultunits.mass)
            @test ed_M.rs isa u.Quantity{Float64}
            ed_Σ0 = ExponentialDisk(1.0*defaultunits.length; Σ0=1.0*defaultunits.surfacedensity)
            @test ed_Σ0.rs isa u.Quantity{Float64}
            # ExponentialDiskDHI constructor
            ed_DHI = ExponentialDiskDHI(1.0*defaultunits.length, 1.0*defaultunits.mass)
            @test ed_DHI.rs isa u.Quantity{Float64}
        end
        @testset "Float32" begin
            ed_plain = ExponentialDisk(1.0f0, 1.0f0)
            @test ed_plain.Σ0 == 1.0f0
            @test ed_plain.rs == 1.0f0
            ed_unitful = ExponentialDisk(1.0f0*defaultunits.surfacedensity, 1.0f0*defaultunits.length)
            @test ed_unitful.Σ0 isa u.Quantity{Float32}
        end
    end

    @testset "Plummer" begin
        @testset "Float64" begin
            # Units are now PRESERVED in struct fields
            plummer_unitful = Plummer(1.0*defaultunits.mass, 1.0*defaultunits.length)
            @test plummer_unitful.M isa u.Quantity{Float64}
            @test plummer_unitful.a isa u.Quantity{Float64}
            # Plain number constructor still works
            plummer_plain = Plummer(1.0, 1.0)
            @test plummer_plain.M == 1.0
            @test plummer_plain.a == 1.0
        end
        @testset "Float32" begin
            plummer_unitful = Plummer(1.0f0*defaultunits.mass, 1.0f0*defaultunits.length)
            @test plummer_unitful.M isa u.Quantity{Float32}
            @test plummer_unitful.a isa u.Quantity{Float32}
        end
    end

    @testset "GeneralIsothermal" begin
        @testset "Float64" begin
            # Units are now PRESERVED in struct fields
            gi_unitful = GeneralIsothermal(1.0*defaultunits.density, 1.0*defaultunits.length, 1.0)
            @test gi_unitful.ρ0 isa u.Quantity{Float64}
            @test gi_unitful.rs isa u.Quantity{Float64}
            @test gi_unitful.α == 1.0  # α is dimensionless
            # Plain number constructor still works
            gi_plain = GeneralIsothermal(1.0, 1.0, 1.0)
            @test gi_plain.ρ0 == 1.0
            @test gi_plain.rs == 1.0
        end
        @testset "Float32" begin
            gi_unitful = GeneralIsothermal(1.0f0*defaultunits.density, 1.0f0*defaultunits.length, 1.0f0)
            @test gi_unitful.ρ0 isa u.Quantity{Float32}
            @test gi_unitful.rs isa u.Quantity{Float32}
        end
    end

    @testset "SIS" begin
        @testset "Float64" begin
            # SIS returns GeneralIsothermal with units PRESERVED
            sis_unitful = SIS(1.0*defaultunits.density, 1.0*defaultunits.length)
            @test sis_unitful isa GeneralIsothermal
            @test sis_unitful.ρ0 isa u.Quantity{Float64}
            @test sis_unitful.rs isa u.Quantity{Float64}
            # Plain number constructor still works  
            sis_plain = SIS(1.0, 1.0)
            @test sis_plain.ρ0 == 1.0
            @test sis_plain.rs == 1.0
            # Test with M, Rmax constructor
            sis_M = SIS(1.0*defaultunits.length, 1.0*defaultunits.mass, 10.0*defaultunits.length)
            @test sis_M.rs isa u.Quantity{Float64}
        end
        @testset "Float32" begin
            sis_unitful = SIS(1.0f0*defaultunits.density, 1.0f0*defaultunits.length)
            @test sis_unitful.ρ0 isa u.Quantity{Float32}
        end
    end    
end


######################################################################
######################################################################


# This testset will exhaustively test the Unitful.jl extension methods. For simplicity we will test with a single density profile (GeneralIsothermal) and a single surface density profile (ExponentialDisk). If the units work on these types, they should work on all others as well.
@testset "GeneralIsothermal Units" begin
    @testset "Float64" begin
        # As tested in the constructor section above,
        # this should construct d = GeneralIsothermal(1.0, 1.0, 1.0)
        d = GeneralIsothermal(1.0 * defaultunits.density,
                              1.0 * defaultunits.length, 1.0)
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
        @test @inferred ρ(defaultunits.density,d,1.0) == 1.0*defaultunits.density
        @test ρ(defaultunits.density,d,1.0) |> u.ustrip isa Float64
        @test ρ(defaultunits.density,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred ρ(defaultunits.density,d,1.0*defaultunits.length) == 1.0*defaultunits.density
        @test ρ(defaultunits.density,d,1.0*defaultunits.length) |> u.ustrip isa Float64
        @test ρ(defaultunits.density,d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        @test @inferred ρ(d,1.0*defaultunits.length) == 1.0*defaultunits.density
        @test ρ(d,1.0*defaultunits.length) |> u.ustrip isa Float64
        @test ρ(d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        ####################################
        @test @inferred invρ(defaultunits.length,d,1.0) == 1.0*defaultunits.length
        @test invρ(defaultunits.length,d,1.0) |> u.ustrip isa Float64
        @test invρ(defaultunits.length,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred invρ(defaultunits.length,d,1.0*defaultunits.density) == 1.0*defaultunits.length
        @test invρ(defaultunits.length,d,1.0*defaultunits.density) |> u.ustrip isa Float64
        @test invρ(defaultunits.length,d,1.0f0*defaultunits.density) |> u.ustrip isa Float64
        @test @inferred invρ(d,1.0*defaultunits.density) == 1.0*defaultunits.length
        @test invρ(d,1.0*defaultunits.density) |> u.ustrip isa Float64
        @test invρ(d,1.0f0*defaultunits.density) |> u.ustrip isa Float64
        ####################################
        @test @inferred ∇ρ(defaultunits.∇density,d,1.0) == -1.0*defaultunits.∇density
        @test ∇ρ(defaultunits.∇density,d,1.0) |> u.ustrip isa Float64
        @test ∇ρ(defaultunits.∇density,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred ∇ρ(defaultunits.∇density,d,1.0*defaultunits.length) == -1.0*defaultunits.∇density
        @test ∇ρ(defaultunits.∇density,d,1.0*defaultunits.length) |> u.ustrip isa Float64
        @test ∇ρ(defaultunits.∇density,d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        @test @inferred ∇ρ(d,1.0*defaultunits.length) == -1.0*defaultunits.∇density
        @test ∇ρ(d,1.0*defaultunits.length) |> u.ustrip isa Float64
        @test ∇ρ(d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        ####################################
        @test @inferred Σ(defaultunits.surfacedensity,d,1.0) == Inf*defaultunits.surfacedensity
        @test Σ(defaultunits.surfacedensity,d,1.0) |> u.ustrip isa Float64
        @test Σ(defaultunits.surfacedensity,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred Σ(defaultunits.surfacedensity,d,1.0*defaultunits.length) == Inf*defaultunits.surfacedensity
        @test Σ(defaultunits.surfacedensity,d,1.0*defaultunits.length) |> u.ustrip isa Float64
        @test Σ(defaultunits.surfacedensity,d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        @test @inferred Σ(d,1.0*defaultunits.length) == Inf*defaultunits.surfacedensity
        @test Σ(d,1.0*defaultunits.length) |> u.ustrip isa Float64        
        @test Σ(d,1.0f0*defaultunits.length) |> u.ustrip isa Float64        

        # This should construct d = GeneralIsothermal(1.0, 1.0, 2.0); then test with Real arguments
        d = GeneralIsothermal(1.0 * defaultunits.density, 1.0 * defaultunits.length, 2.0)
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
        @test @inferred Mproj(d,1.0) ≈ 2*π^2
        @test Mproj(d,1.0) isa Float64
        @test @inferred Mproj(d,1.0f0) ≈ 2*π^2
        @test Mproj(d,1.0f0) isa Float64
        @test @inferred ∇Mproj(d,1.0) ≈ 2*π^2
        @test ∇Mproj(d,1.0) isa Float64
        @test @inferred ∇Mproj(d,1.0f0) ≈ 2*π^2
        @test ∇Mproj(d,1.0f0) isa Float64
        @test @inferred invMproj(d,2*π^2) ≈ 1.0
        @test invMproj(d,2*π^2) isa Float64
        @test @inferred invMproj(d,Float32(2*π^2)) ≈ 1.0f0
        @test invMproj(d,Float32(2*π^2)) isa Float64
        @test @inferred dynamical_time(d,3.0) ≈ 6.267613877418796e11
        @test dynamical_time(d,3.0) isa Float64
        @test @inferred dynamical_time(d,3.0f0) ≈ 6.267613877418796f11
        @test dynamical_time(d,3.0f0) isa Float64
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
        @test @inferred invΣ(defaultunits.length,d,π) == 1.0*defaultunits.length
        @test invΣ(defaultunits.length,d,π) |> u.ustrip isa Float64
        @test invΣ(defaultunits.length,d,Float32(π)) |> u.ustrip isa Float64
        @test @inferred invΣ(defaultunits.length,d,π*defaultunits.surfacedensity) == 1*defaultunits.length
        @test invΣ(defaultunits.length,d,π*defaultunits.surfacedensity) |> u.ustrip isa Float64
        @test invΣ(defaultunits.length,d,Float32(π)*defaultunits.surfacedensity) |> u.ustrip isa Float64
        @test @inferred invΣ(d,π*defaultunits.surfacedensity) == 1.0*defaultunits.length
        @test invΣ(d,π*defaultunits.surfacedensity) |> u.ustrip isa Float64 
        @test invΣ(d,Float32(π)*defaultunits.surfacedensity) |> u.ustrip isa Float64
        ####################################
        @test @inferred ∇Σ(defaultunits.density,d,1.0) == -π*defaultunits.density
        @test ∇Σ(defaultunits.density,d,1.0) |> u.ustrip isa Float64
        @test ∇Σ(defaultunits.density,d,Float32(1.0)) |> u.ustrip isa Float64
        @test @inferred ∇Σ(defaultunits.density,d,1.0*defaultunits.length) == -π*defaultunits.density
        @test ∇Σ(defaultunits.density,d,1.0*defaultunits.length) |> u.ustrip isa Float64
        @test ∇Σ(defaultunits.density,d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        @test @inferred ∇Σ(d,1.0*defaultunits.length) == -π*defaultunits.density
        @test ∇Σ(d,1.0*defaultunits.length) |> u.ustrip isa Float64 
        @test ∇Σ(d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        ####################################
        @test @inferred M(defaultunits.mass,d,1.0) == 4π*defaultunits.mass
        @test M(defaultunits.mass,d,1.0) |> u.ustrip isa Float64
        @test M(defaultunits.mass,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred M(defaultunits.mass,d,1.0*defaultunits.length) == 4π*defaultunits.mass
        @test M(defaultunits.mass,d,1.0*defaultunits.length) |> u.ustrip isa Float64
        @test M(defaultunits.mass,d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        @test @inferred M(d,1.0*defaultunits.length) == 4π*defaultunits.mass
        @test M(d,1.0*defaultunits.length) |> u.ustrip isa Float64
        @test M(d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        ####################################
        @test @inferred invM(defaultunits.length,d,4π) == 1.0*defaultunits.length
        @test invM(defaultunits.length,d,4π) |> u.ustrip isa Float64
        @test invM(defaultunits.length,d,Float32(4π)) |> u.ustrip isa Float64
        @test @inferred invM(defaultunits.length,d,4π*defaultunits.mass) == 1.0*defaultunits.length
        @test invM(defaultunits.length,d,4π*defaultunits.mass) |> u.ustrip isa Float64
        @test invM(defaultunits.length,d,Float32(4π)*defaultunits.mass) |> u.ustrip isa Float64
        @test @inferred invM(d,4π*defaultunits.mass) == 1.0*defaultunits.length
        @test invM(d,4π*defaultunits.mass) |> u.ustrip isa Float64
        @test invM(d,Float32(4π)*defaultunits.mass) |> u.ustrip isa Float64
        ####################################
        @test @inferred ∇M(defaultunits.∇mass,d,1.0) == 4π*defaultunits.∇mass
        @test ∇M(defaultunits.∇mass,d,1.0) |> u.ustrip isa Float64
        @test ∇M(defaultunits.∇mass,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred ∇M(defaultunits.∇mass,d,1.0*defaultunits.length) == 4π*defaultunits.∇mass
        @test ∇M(defaultunits.∇mass,d,1.0*defaultunits.length) |> u.ustrip isa Float64
        @test ∇M(defaultunits.∇mass,d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        @test @inferred ∇M(d,1.0*defaultunits.length) == 4π*defaultunits.∇mass
        @test ∇M(d,1.0*defaultunits.length) |> u.ustrip isa Float64
        @test ∇M(d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        ####################################
        @test @inferred Mproj(defaultunits.mass,d,1.0) ≈ 2*π^2 * defaultunits.mass
        @test Mproj(defaultunits.mass,d,1.0) |> u.ustrip isa Float64
        @test Mproj(defaultunits.mass,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred Mproj(defaultunits.mass,d,1.0*defaultunits.length) ≈ 2*π^2 * defaultunits.mass
        @test Mproj(defaultunits.mass,d,1.0*defaultunits.length) |> u.ustrip isa Float64
        @test Mproj(defaultunits.mass,d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        @test @inferred Mproj(d,1.0*defaultunits.length) ≈ 2*π^2 * defaultunits.mass
        @test Mproj(d,1.0*defaultunits.length) |> u.ustrip isa Float64
        @test Mproj(d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        ####################################
        @test @inferred ∇Mproj(defaultunits.∇mass,d,1.0) ≈ 2*π^2*defaultunits.∇mass
        @test ∇Mproj(defaultunits.∇mass,d,1.0) |> u.ustrip isa Float64
        @test ∇Mproj(defaultunits.∇mass,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred ∇Mproj(defaultunits.∇mass,d,1.0*defaultunits.length) ≈ 2*π^2*defaultunits.∇mass
        @test ∇Mproj(defaultunits.∇mass,d,1.0*defaultunits.length) |> u.ustrip isa Float64
        @test ∇Mproj(defaultunits.∇mass,d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        @test @inferred ∇Mproj(d,1.0*defaultunits.length) ≈ 2*π^2*defaultunits.∇mass
        @test ∇Mproj(d,1.0*defaultunits.length) |> u.ustrip isa Float64
        @test ∇Mproj(d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        ####################################
        @test @inferred invMproj(defaultunits.length,d,2*π^2) ≈ 1.0*defaultunits.length
        @test invMproj(defaultunits.length,d,2*π^2) |> u.ustrip isa Float64
        @test invMproj(defaultunits.length,d,Float32(2*π^2)) |> u.ustrip isa Float64
        @test @inferred invMproj(defaultunits.length,d,2*π^2*defaultunits.mass) ≈ 1.0*defaultunits.length
        @test invMproj(defaultunits.length,d,2*π^2*defaultunits.mass) |> u.ustrip isa Float64
        @test invMproj(defaultunits.length,d,Float32(2*π^2)*defaultunits.mass) |> u.ustrip isa Float64
        @test @inferred invMproj(d,2*π^2*defaultunits.mass) ≈ 1.0*defaultunits.length
        @test invMproj(d,2*π^2*defaultunits.mass) |> u.ustrip isa Float64
        @test invMproj(d,Float32(2*π^2)*defaultunits.mass) |> u.ustrip isa Float64
        ####################################
        @test @inferred dynamical_time(defaultunits.time,d,3.0) ≈ 6.267613877418796e11*defaultunits.time
        @test dynamical_time(defaultunits.time,d,3.0) |> u.ustrip isa Float64
        @test dynamical_time(defaultunits.time,d,Float32(3.0)) |> u.ustrip isa Float64
        @test @inferred dynamical_time(defaultunits.time,d,3.0*defaultunits.length) ≈ 6.267613877418796e11*defaultunits.time
        @test dynamical_time(defaultunits.time,d,3.0*defaultunits.length) |> u.ustrip isa Float64
        @test dynamical_time(defaultunits.time,d,Float32(3.0)*defaultunits.length) |> u.ustrip isa Float64
        @test @inferred dynamical_time(d,3.0*defaultunits.length) ≈ 6.267613877418796e11*defaultunits.time
        @test dynamical_time(d,3.0*defaultunits.length) |> u.ustrip isa Float64
        @test dynamical_time(d,Float32(3.0)*defaultunits.length) |> u.ustrip isa Float64
        ####################################
        @test @inferred Vcirc(defaultunits.velocity,d,1.0) == Vcirc(d,1.0) * defaultunits.velocity
        @test Vcirc(defaultunits.velocity,d,1.0) |> u.ustrip isa Float64
        @test Vcirc(defaultunits.velocity,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred Vcirc(defaultunits.velocity,d,1.0*defaultunits.length) == Vcirc(d,1.0) * defaultunits.velocity
        @test Vcirc(defaultunits.velocity,d,1.0*defaultunits.length) |> u.ustrip isa Float64
        @test Vcirc(defaultunits.velocity,d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        @test @inferred Vcirc(d,1.0*defaultunits.length) == Vcirc(d,1.0) * defaultunits.velocity
        @test Vcirc(d,1.0*defaultunits.length) |> u.ustrip isa Float64
        @test Vcirc(d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        ####################################
        @test_throws DomainError Vesc(defaultunits.velocity,d,1.0)   # α=2 is not valid for Vesc
        @test_throws DomainError Vesc(defaultunits.velocity,d,1.0*defaultunits.length)   # α=2 is not valid for Vesc
        @test_throws DomainError Vesc(d,1.0*defaultunits.length)   # α=2 is not valid for Vesc
        ####################################
        let d = GeneralIsothermal(1.0*defaultunits.density,1.0*defaultunits.length,2.5) # introduce new d just for Vesc
            @test d == GeneralIsothermal(1.0,1.0,2.5)
            @test @inferred Vesc(defaultunits.velocity,d,1.0) == Vesc(d,1.0) * defaultunits.velocity
            @test Vesc(defaultunits.velocity,d,1.0) |> u.ustrip isa Float64
            @test Vesc(defaultunits.velocity,d,1.0f0) |> u.ustrip isa Float64
            @test @inferred Vesc(defaultunits.velocity,d,1.0*defaultunits.length) == Vesc(d,1.0) * defaultunits.velocity
            @test Vesc(defaultunits.velocity,d,1.0*defaultunits.length) |> u.ustrip isa Float64
            @test Vesc(defaultunits.velocity,d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
            @test @inferred Vesc(d,1.0*defaultunits.length) == Vesc(d,1.0) * defaultunits.velocity
            @test Vesc(d,1.0*defaultunits.length) |> u.ustrip isa Float64
            @test Vesc(d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        end
        ####################################
        # begin # potentials; for some reason these mess up Emacs julia mode so stick them in a begin block
        @test @inferred Φ(defaultunits.Φunit,d,1.0) == Φ(d,1.0) * defaultunits.Φunit
        @test Φ(defaultunits.Φunit,d,1.0) |> u.ustrip isa Float64
        @test Φ(defaultunits.Φunit,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred Φ(defaultunits.Φunit,d,1.0*defaultunits.length) == Φ(d,1.0) * defaultunits.Φunit
        @test Φ(defaultunits.Φunit,d,1.0*defaultunits.length) |> u.ustrip isa Float64
        @test Φ(defaultunits.Φunit,d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        @test @inferred Φ(d,1.0*defaultunits.length) == Φ(d,1.0) * defaultunits.Φunit
        @test Φ(d,1.0*defaultunits.length) |> u.ustrip isa Float64
        @test Φ(d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        ####################################
        @test @inferred ∇Φ(defaultunits.∇Φunit,d,1.0) == ∇Φ(d,1.0) * defaultunits.∇Φunit
        @test @inferred ∇Φ(defaultunits.∇Φunit,d,1.0) == ∇Φ(d,1.0) * defaultunits.∇Φunit
        @test ∇Φ(defaultunits.∇Φunit,d,1.0) |> u.ustrip isa Float64
        @test ∇Φ(defaultunits.∇Φunit,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred ∇Φ(defaultunits.∇Φunit,d,1.0*defaultunits.length) == ∇Φ(d,1.0) * defaultunits.∇Φunit
        @test ∇Φ(defaultunits.∇Φunit,d,1.0*defaultunits.length) |> u.ustrip isa Float64
        @test ∇Φ(defaultunits.∇Φunit,d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        @test @inferred ∇Φ(d,1.0*defaultunits.length) == ∇Φ(d,1.0) * defaultunits.∇Φunit
        @test ∇Φ(d,1.0*defaultunits.length) |> u.ustrip isa Float64
        @test ∇Φ(d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        ####################################
        @test @inferred ∇∇Φ(defaultunits.∇∇Φunit,d,1.0) == ∇∇Φ(d,1.0) * defaultunits.∇∇Φunit
        @test ∇∇Φ(defaultunits.∇∇Φunit,d,1.0) |> u.ustrip isa Float64
        @test ∇∇Φ(defaultunits.∇∇Φunit,d,1.0f0) |> u.ustrip isa Float64
        @test @inferred ∇∇Φ(defaultunits.∇∇Φunit,d,1.0*defaultunits.length) == ∇∇Φ(d,1.0) * defaultunits.∇∇Φunit
        @test ∇∇Φ(defaultunits.∇∇Φunit,d,1.0*defaultunits.length) |> u.ustrip isa Float64
        @test ∇∇Φ(defaultunits.∇∇Φunit,d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        @test @inferred ∇∇Φ(d,1.0*defaultunits.length) == ∇∇Φ(d,1.0) * defaultunits.∇∇Φunit
        @test ∇∇Φ(d,1.0*defaultunits.length) |> u.ustrip isa Float64
        @test ∇∇Φ(d,1.0f0*defaultunits.length) |> u.ustrip isa Float64
        ####################################
    end
end
