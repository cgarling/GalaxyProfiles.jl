# New tests for Unitful quantity preservation
using GalaxyProfiles
using Test
import Unitful as u
import UnitfulAstro as ua

@testset "Unitful Quantity Preservation" begin
    
    @testset "Plummer with units preserved" begin
        # Create Plummer with Unitful quantities
        d = Plummer(1e8*ua.Msun, 0.5*ua.kpc)
        
        # Check that fields preserve units
        @test d.M isa u.Quantity
        @test d.a isa u.Quantity
        @test u.dimension(d.M) == u.dimension(1.0*ua.Msun)
        @test u.dimension(d.a) == u.dimension(1.0*ua.kpc)
        
        # Test that methods return correct units
        @test ρ(d, 0.5*ua.kpc) isa u.Quantity
        @test u.dimension(ρ(d, 0.5*ua.kpc)) == u.dimension(1.0*ua.Msun/ua.kpc^3)
        
        @test Σ(d, 0.5*ua.kpc) isa u.Quantity
        @test u.dimension(Σ(d, 0.5*ua.kpc)) == u.dimension(1.0*ua.Msun/ua.kpc^2)
        
        @test M(d, 1.0*ua.kpc) isa u.Quantity
        @test u.dimension(M(d, 1.0*ua.kpc)) == u.dimension(1.0*ua.Msun)
        
        @test Mproj(d, 1.0*ua.kpc) isa u.Quantity
        @test u.dimension(Mproj(d, 1.0*ua.kpc)) == u.dimension(1.0*ua.Msun)
        
        # Test parameter accessors
        @test Mtot(d) == 1e8*ua.Msun
        @test scale_radius(d) == 0.5*ua.kpc
        
        # TODO: Test dynamical methods when they're updated
        # @test Φ(d, 1.0*ua.kpc) isa u.Quantity
        # @test ∇Φ(d, 1.0*ua.kpc) isa u.Quantity
        # @test ∇∇Φ(d, 1.0*ua.kpc) isa u.Quantity
    end
    
    @testset "NFW with units preserved" begin
        # Create NFW with Unitful quantities
        d = NFW(1e6*ua.Msun/ua.kpc^3, 10.0*ua.kpc)
        
        # Check that fields preserve units
        @test d.ρ0 isa u.Quantity
        @test d.rs isa u.Quantity
        @test u.dimension(d.ρ0) == u.dimension(1.0*ua.Msun/ua.kpc^3)
        @test u.dimension(d.rs) == u.dimension(1.0*ua.kpc)
        
        # Test that methods return correct units
        @test ρ(d, 5.0*ua.kpc) isa u.Quantity
        @test u.dimension(ρ(d, 5.0*ua.kpc)) == u.dimension(1.0*ua.Msun/ua.kpc^3)
        
        @test M(d, 10.0*ua.kpc) isa u.Quantity
        @test u.dimension(M(d, 10.0*ua.kpc)) == u.dimension(1.0*ua.Msun)
    end
    
    @testset "GeneralIsothermal with units preserved" begin
        # Create GeneralIsothermal with Unitful quantities
        d = GeneralIsothermal(1e6*ua.Msun/ua.kpc^3, 5.0*ua.kpc, 2.0)
        
        # Check that fields preserve units
        @test d.ρ0 isa u.Quantity
        @test d.rs isa u.Quantity
        @test d.α isa Real  # α is dimensionless
        @test u.dimension(d.ρ0) == u.dimension(1.0*ua.Msun/ua.kpc^3)
        @test u.dimension(d.rs) == u.dimension(1.0*ua.kpc)
        
        # Test that methods return correct units
        @test ρ(d, 3.0*ua.kpc) isa u.Quantity
        @test u.dimension(ρ(d, 3.0*ua.kpc)) == u.dimension(1.0*ua.Msun/ua.kpc^3)
    end
    
    @testset "Sampling functions with units" begin
        # Create Plummer with Unitful quantities
        d = Plummer(1e8*ua.Msun, 1.0*ua.kpc)
        
        # TODO: Sampling functions need to be updated to handle Unitful quantities
        # The current implementation tries to assign Unitful quantities to Float64 arrays
        
        # For now, test that quantile functions work
        @test quantile2D(d, 0.5) isa u.Quantity
        @test u.dimension(quantile2D(d, 0.5)) == u.dimension(1.0*ua.kpc)
        
        @test quantile3D(d, 0.5) isa u.Quantity
        @test u.dimension(quantile3D(d, 0.5)) == u.dimension(1.0*ua.kpc)
        
        # Test that CDFs work correctly
        @test cdf2D(d, quantile2D(d, 0.5)) ≈ 0.5
        @test cdf3D(d, quantile3D(d, 0.5)) ≈ 0.5
    end
    
    @testset "Mixed units input" begin
        # Test with different length units
        d1 = Plummer(1e8*ua.Msun, 1.0*ua.kpc)
        d2 = Plummer(1e8*ua.Msun, 1000.0*ua.pc)
        
        # These should give same results when we compute ρ at the same physical radius
        ρ1 = ρ(d1, 0.5*ua.kpc)
        ρ2 = ρ(d2, 500.0*ua.pc)
        @test u.uconvert(ua.Msun/ua.kpc^3, ρ1) ≈ u.uconvert(ua.Msun/ua.kpc^3, ρ2)
    end
    
    @testset "Backwards compatibility with plain numbers" begin
        # Test that plain numbers still work
        d = Plummer(1.0, 1.0)
        @test d.M == 1.0
        @test d.a == 1.0
        @test d.M isa Float64
        @test d.a isa Float64
        
        # Methods should return plain numbers
        @test ρ(d, 1.0) isa Float64
        @test Σ(d, 1.0) isa Float64
        @test M(d, 1.0) isa Float64
        
        # Sampling should return plain numbers
        r_samples = sample2D_r(d, 10)
        @test all(r isa Float64 for r in r_samples)
    end
    
end
