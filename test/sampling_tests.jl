using GalaxyProfiles
using Test
using Random: MersenneTwister
using Statistics: mean, median
using SkyCoords: ICRSCoords, lon, lat

@testset "Random Sampling Methods" begin
    
    # Set up test profile
    plummer = Plummer(1000.0, 0.5)
    
    # Set up test coordinates and distance for 2D sampling
    center = ICRSCoords(deg2rad(270.0), deg2rad(66.0))
    distance = 100.0  # kpc
    
    # Create a deterministic RNG for reproducible tests
    rng = MersenneTwister(42)
    
    @testset "sample2D_r!" begin
        @testset "Basic functionality" begin
            # Test with pre-allocated array
            x = zeros(100)
            result = sample2D_r!(rng, plummer, x)
            
            # Should mutate the input array
            @test result === x
            @test all(x .> 0)  # All radii should be positive
            @test all(isfinite.(x))  # All values should be finite
            
            # Test without explicit rng (uses default_rng internally)
            y = zeros(50)
            result2 = sample2D_r!(plummer, y)
            @test result2 === y
            @test all(y .> 0)
        end
        
        @testset "Different array shapes" begin
            # Test with different array dimensions
            x_1d = zeros(10)
            sample2D_r!(rng, plummer, x_1d)
            @test length(x_1d) == 10
            @test all(x_1d .> 0)
            
            x_2d = zeros(5, 5)
            sample2D_r!(rng, plummer, x_2d)
            @test size(x_2d) == (5, 5)
            @test all(x_2d .> 0)
            
            x_3d = zeros(2, 3, 4)
            sample2D_r!(rng, plummer, x_3d)
            @test size(x_3d) == (2, 3, 4)
            @test all(x_3d .> 0)
        end
        
        @testset "Type stability" begin
            x_f64 = zeros(Float64, 10)
            sample2D_r!(rng, plummer, x_f64)
            @test eltype(x_f64) == Float64
            
            x_f32 = zeros(Float32, 10)
            sample2D_r!(rng, plummer, x_f32)
            @test eltype(x_f32) == Float32
        end
    end
    
    @testset "sample2D_r" begin
        @testset "Single sample" begin
            # Test single sample
            r1 = sample2D_r(rng, plummer)
            @test r1 isa Real
            @test r1 > 0
            @test isfinite(r1)
            
            # Test without explicit rng
            r2 = sample2D_r(plummer)
            @test r2 isa Real
            @test r2 > 0
        end
        
        @testset "Multiple samples with Vararg" begin
            # Test with single integer
            samples_1d = sample2D_r(rng, plummer, 100)
            @test samples_1d isa Vector
            @test length(samples_1d) == 100
            @test all(samples_1d .> 0)
            
            # Test with multiple integers (2D array)
            samples_2d = sample2D_r(rng, plummer, 10, 20)
            @test samples_2d isa Matrix
            @test size(samples_2d) == (10, 20)
            @test all(samples_2d .> 0)
            
            # Test with three integers (3D array)
            samples_3d = sample2D_r(rng, plummer, 2, 3, 4)
            @test samples_3d isa Array{<:Real, 3}
            @test size(samples_3d) == (2, 3, 4)
            @test all(samples_3d .> 0)
        end
        
        @testset "Multiple samples with Dims" begin
            # Test with Dims (tuple)
            samples_tuple = sample2D_r(rng, plummer, (10, 20))
            @test samples_tuple isa Matrix
            @test size(samples_tuple) == (10, 20)
            @test all(samples_tuple .> 0)
            
            # Test without explicit rng
            samples_norng = sample2D_r(plummer, (5, 5))
            @test size(samples_norng) == (5, 5)
        end
        
        @testset "Distribution properties" begin
            # Implement here some check that the samples follow the expected distribution
        end
    end
    
    @testset "sample2D (sky coordinates)" begin
        @testset "Single sample" begin
            # Test single sample
            coord = sample2D(rng, plummer, center, distance)
            @test coord isa ICRSCoords
            
            # Coordinates should be close to center for small radii
            # (though with some scatter)
            @test isfinite(lon(coord))
            @test isfinite(lat(coord))
            
            # Test without explicit rng
            coord2 = sample2D(plummer, center, distance)
            @test coord2 isa ICRSCoords
        end
        
        @testset "Multiple samples with Vararg" begin
            # Test with single integer
            coords_1d = sample2D(rng, plummer, center, distance, 50)
            @test coords_1d isa Vector{<:ICRSCoords}
            @test length(coords_1d) == 50
            
            # Test with multiple integers
            coords_2d = sample2D(rng, plummer, center, distance, 5, 10)
            @test coords_2d isa Matrix{<:ICRSCoords}
            @test size(coords_2d) == (5, 10)
        end
        
        @testset "Multiple samples with Dims" begin
            # Test with Dims (tuple)
            coords_tuple = sample2D(rng, plummer, center, distance, (5, 10))
            @test coords_tuple isa Matrix{<:ICRSCoords}
            @test size(coords_tuple) == (5, 10)
            
            # Test without explicit rng
            coords_norng = sample2D(plummer, center, distance, (3, 3))
            @test coords_norng isa Matrix{<:ICRSCoords}
            @test size(coords_norng) == (3, 3)
        end
        
        @testset "Angular separation calculation" begin
            # Some validation that the angular separations correspond to the sampled radii
        end
        
        @testset "Position angle distribution" begin
            # Sample many points and check that position angles are uniform
            coords = sample2D(rng, plummer, center, distance, 1000)
            
            # Calculate position angles relative to center
            pas = [begin
                Δlon = lon(c) - lon(center)
                Δlat = lat(c) - lat(center)
                atan(Δlon, Δlat)
            end for c in coords]
            
            # Bin the position angles and check for uniformity
            # This is a rough test - with 1000 samples, each quadrant should have ~250
            n_q1 = count(0 <= pa < π/2 for pa in pas)
            n_q2 = count(π/2 <= pa < π for pa in pas)
            n_q3 = count(-π <= pa < -π/2 for pa in pas)
            n_q4 = count(-π/2 <= pa < 0 for pa in pas)
            
            # Each quadrant should have roughly 250 ± sqrt(250) ≈ 250 ± 16
            # Using 3σ bounds: roughly 200-300 per quadrant
            @test 150 < n_q1 < 350
            @test 150 < n_q2 < 350
            @test 150 < n_q3 < 350
            @test 150 < n_q4 < 350
        end
    end
    
    @testset "sample3D_r!" begin
        @testset "Basic functionality" begin
            # Test with pre-allocated array
            x = zeros(100)
            result = sample3D_r!(rng, plummer, x)
            
            # Should mutate the input array
            @test result === x
            @test all(x .> 0)  # All radii should be positive
            @test all(isfinite.(x))  # All values should be finite
            
            # Test without explicit rng
            y = zeros(50)
            result2 = sample3D_r!(plummer, y)
            @test result2 === y
            @test all(y .> 0)
        end
        
        @testset "Different array shapes" begin
            # Test with different array dimensions
            x_1d = zeros(10)
            sample3D_r!(rng, plummer, x_1d)
            @test length(x_1d) == 10
            @test all(x_1d .> 0)
            
            x_2d = zeros(5, 5)
            sample3D_r!(rng, plummer, x_2d)
            @test size(x_2d) == (5, 5)
            @test all(x_2d .> 0)
        end
    end
    
    @testset "sample3D_r" begin
        @testset "Single sample" begin
            # Test single sample
            r1 = sample3D_r(rng, plummer)
            @test r1 isa Real
            @test r1 > 0
            @test isfinite(r1)
            
            # Test without explicit rng
            r2 = sample3D_r(plummer)
            @test r2 isa Real
            @test r2 > 0
        end
        
        @testset "Multiple samples with Vararg" begin
            # Test with single integer
            samples_1d = sample3D_r(rng, plummer, 100)
            @test samples_1d isa Vector
            @test length(samples_1d) == 100
            @test all(samples_1d .> 0)
            
            # Test with multiple integers (2D array)
            samples_2d = sample3D_r(rng, plummer, 10, 20)
            @test samples_2d isa Matrix
            @test size(samples_2d) == (10, 20)
            @test all(samples_2d .> 0)
        end
        
        @testset "Multiple samples with Dims" begin
            # Test with Dims (tuple)
            samples_tuple = sample3D_r(rng, plummer, (10, 20))
            @test samples_tuple isa Matrix
            @test size(samples_tuple) == (10, 20)
            @test all(samples_tuple .> 0)
            
            # Test without explicit rng
            samples_norng = sample3D_r(plummer, (5, 5))
            @test size(samples_norng) == (5, 5)
        end
        
        @testset "Distribution properties" begin
            # Some check that the 3D samples follow the expected distribution
        end
        
        @testset "3D vs 2D comparison" begin
            # 3D samples should be more extended than 2D samples on average
            n_samples = 5000
            samples_2d = sample2D_r(rng, plummer, n_samples)
            samples_3d = sample3D_r(rng, plummer, n_samples)
            
            # Median 3D radius should be larger
            @test median(samples_3d) > median(samples_2d)
            
            # Mean 3D radius should be larger
            @test mean(samples_3d) > mean(samples_2d)
        end
    end
    
    @testset "Edge cases" begin
        @testset "Small number of samples" begin
            # Should work with just 1 sample
            @test sample2D_r(rng, plummer, 1) isa Vector
            @test length(sample2D_r(rng, plummer, 1)) == 1
            
            @test sample3D_r(rng, plummer, 1) isa Vector
            @test length(sample3D_r(rng, plummer, 1)) == 1
        end
        
        @testset "Different distances for sample2D" begin
            # Very small distance
            coords_near = sample2D(rng, plummer, center, 1.0, 10)
            @test all(c isa ICRSCoords for c in coords_near)
            
            # Very large distance
            coords_far = sample2D(rng, plummer, center, 1e9, 10)
            @test all(c isa ICRSCoords for c in coords_far)
        end
    end
end
