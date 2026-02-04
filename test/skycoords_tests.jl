using GalaxyProfiles
using Test
using SkyCoords
using Random

# Import the function from the extension
if isdefined(Base, :get_extension)
    ext = Base.get_extension(GalaxyProfiles, :GalaxyProfilesSkyCoords)
    const sample2D_skycoords = ext.sample2D_skycoords
else
    # For Julia < 1.9
    const sample2D_skycoords = GalaxyProfiles.GalaxyProfilesSkyCoords.sample2D_skycoords
end

@testset "SkyCoords Extension" begin
    
    @testset "Basic functionality" begin
        # Create a Plummer profile (has quantile2D implemented)
        d = Plummer(1000.0, 10.0)  # M=1000, a=10.0 kpc
        
        # Create galaxy center at RA=0, Dec=0 (in radians)
        center = ICRSCoords(0.0, 0.0)
        
        # Distance to galaxy in kpc (matching profile units)
        distance = 10000.0  # 10 Mpc = 10000 kpc
        
        # Test single sample
        rng = MersenneTwister(42)
        coord = sample2D_skycoords(rng, d, center, distance)
        @test coord isa ICRSCoords
        @test coord.ra isa Float64
        @test coord.dec isa Float64
        
        # Test that coordinates are reasonable (not too far from center)
        # Maximum reasonable offset would be a few times the scale radius
        max_reasonable_sep = 5 * 10.0 / distance  # 5 * scale_radius / distance in radians
        @test abs(coord.ra) < max_reasonable_sep
        @test abs(coord.dec) < max_reasonable_sep
    end
    
    @testset "Multiple samples" begin
        d = Plummer(1000.0, 10.0)  # M=1000, a=10 kpc
        center = ICRSCoords(0.0, 0.0)
        distance = 10000.0
        
        # Test array of samples
        rng = MersenneTwister(42)
        coords = sample2D_skycoords(rng, d, center, distance, 10)
        @test coords isa Vector
        @test length(coords) == 10
        @test all(c -> c isa ICRSCoords, coords)
        
        # Test 2D array
        rng = MersenneTwister(42)
        coords_2d = sample2D_skycoords(rng, d, center, distance, 5, 4)
        @test coords_2d isa Matrix
        @test size(coords_2d) == (5, 4)
        @test all(c -> c isa ICRSCoords, coords_2d)
        
        # Test with Dims tuple
        rng = MersenneTwister(42)
        coords_dims = sample2D_skycoords(rng, d, center, distance, (3, 3))
        @test coords_dims isa Matrix
        @test size(coords_dims) == (3, 3)
        @test all(c -> c isa ICRSCoords, coords_dims)
    end
    
    @testset "No explicit RNG" begin
        d = Plummer(1000.0, 10.0)
        center = ICRSCoords(0.0, 0.0)
        distance = 10000.0
        
        # Single sample without RNG
        coord = sample2D_skycoords(d, center, distance)
        @test coord isa ICRSCoords
        
        # Multiple samples without RNG
        coords = sample2D_skycoords(d, center, distance, 5)
        @test coords isa Vector
        @test length(coords) == 5
    end
    
    @testset "Different coordinate types" begin
        d = Plummer(1000.0, 10.0)
        distance = 10000.0
        rng = MersenneTwister(42)
        
        # Test with GalCoords
        center_gal = GalCoords(0.0, 0.0)
        coord_gal = sample2D_skycoords(rng, d, center_gal, distance)
        @test coord_gal isa GalCoords
    end
    
    @testset "Different profile types" begin
        center = ICRSCoords(0.0, 0.0)
        distance = 10000.0
        rng = MersenneTwister(42)
        
        # Test with Plummer profile
        p = Plummer(1000.0, 5.0)  # M=1000, a=5 kpc
        coord_p = sample2D_skycoords(rng, p, center, distance)
        @test coord_p isa ICRSCoords
        
        # Test with ExponentialDisk (has quantile2D implemented)
        disk = ExponentialDisk(1.0, 3.0)  # Σ0=1.0, rs=3.0 kpc
        coord_disk = sample2D_skycoords(rng, disk, center, distance)
        @test coord_disk isa ICRSCoords
    end
    
    @testset "Reasonable angular separations" begin
        # Create profile and sample
        d = Plummer(1000.0, 10.0)  # M=1000, a = 10 kpc
        center = ICRSCoords(0.0, 0.0)
        distance = 10000.0  # 10 Mpc
        
        rng = MersenneTwister(42)
        coords = sample2D_skycoords(rng, d, center, distance, 100)
        
        # Calculate separations from center
        separations = [SkyCoords.separation(center, c) for c in coords]
        
        # Separations should be positive
        @test all(s -> s >= 0, separations)
        
        # Most samples should be within reasonable distance
        # For Plummer, the distribution has a long tail, so use a large multiple
        max_expected_sep = 50 * 10.0 / distance  # 50 * a / distance in radians (very generous)
        @test all(s -> s < max_expected_sep, separations)
        
        # At least some samples should be non-zero
        @test any(s -> s > 0, separations)
    end
    
    @testset "Consistency check" begin
        # Sample twice with same RNG seed should give same results
        d = Plummer(1000.0, 10.0)
        center = ICRSCoords(0.0, 0.0)
        distance = 10000.0
        
        rng1 = MersenneTwister(12345)
        coords1 = sample2D_skycoords(rng1, d, center, distance, 5)
        
        rng2 = MersenneTwister(12345)
        coords2 = sample2D_skycoords(rng2, d, center, distance, 5)
        
        # Should get identical results
        for (c1, c2) in zip(coords1, coords2)
            @test c1.ra ≈ c2.ra
            @test c1.dec ≈ c2.dec
        end
    end
    
    @testset "Non-zero center" begin
        # Test with non-zero center position
        d = Plummer(1000.0, 10.0)
        # Center at RA=3h (45°), Dec=30°
        center = ICRSCoords(deg2rad(45.0), deg2rad(30.0))
        distance = 10000.0
        
        rng = MersenneTwister(42)
        coords = sample2D_skycoords(rng, d, center, distance, 10)
        
        @test all(c -> c isa ICRSCoords, coords)
        
        # All coords should be near the center
        for coord in coords
            sep = SkyCoords.separation(center, coord)
            # Reasonable separation (within ~10 scale radii)
            @test sep < 10 * 10.0 / distance
        end
    end
end
