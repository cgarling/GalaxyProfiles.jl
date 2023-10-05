# Just test constructors. Due to how we are overloading methods in units.jl via metaprogramming, if the overloaded methods work for one type they should work on others, as long as the call signature is the same. We are therefore not going to retest them all here, as they are mostly covered in general_isothermal_units.jl.
@testset "ExponentialDisk Units" begin
    @testset "Float64" begin
        @test ExponentialDisk(1.0*GalaxyProfiles.defaultunits.surfacedensity,1.0*GalaxyProfiles.defaultunits.length) isa ExponentialDisk{Float64}
        @test ExponentialDisk(1.0*GalaxyProfiles.defaultunits.surfacedensity,1.0*GalaxyProfiles.defaultunits.length) == ExponentialDisk(1.0,1.0)
        @test ExponentialDisk(1.0*GalaxyProfiles.defaultunits.length; M=1.0*GalaxyProfiles.defaultunits.mass) isa ExponentialDisk{Float64}
        @test ExponentialDisk(1.0*GalaxyProfiles.defaultunits.length; Σ0=1.0*GalaxyProfiles.defaultunits.surfacedensity) isa ExponentialDisk{Float64}
        @test ExponentialDiskDHI(1.0*GalaxyProfiles.defaultunits.length,1.0*GalaxyProfiles.defaultunits.mass) isa ExponentialDisk{Float64}
        #################################
        @test ExponentialDiskDHI(1.0*GalaxyProfiles.defaultunits.length,1.0*GalaxyProfiles.defaultunits.mass,1e6*GalaxyProfiles.defaultunits.surfacedensity) isa ExponentialDisk{Float64}
        @test ExponentialDiskDHI(1.0*GalaxyProfiles.defaultunits.length,1.0f0*GalaxyProfiles.defaultunits.mass,1e6*GalaxyProfiles.defaultunits.surfacedensity) isa ExponentialDisk{Float64}        
    end
    @testset "Float32" begin
        @test ExponentialDisk(1.0f0,1.0f0) isa ExponentialDisk{Float32}
        @test ExponentialDisk(1.0f0,1) isa ExponentialDisk{Float32}
        @test ExponentialDisk(1,1.0f0) isa ExponentialDisk{Float32}
        @test ExponentialDisk(1.0f0*GalaxyProfiles.defaultunits.length; M=1.0f0*GalaxyProfiles.defaultunits.mass) isa ExponentialDisk{Float32}
        @test ExponentialDisk(1.0f0*GalaxyProfiles.defaultunits.length; Σ0=1.0f0*GalaxyProfiles.defaultunits.surfacedensity) isa ExponentialDisk{Float32}
        ##################################
        @test ExponentialDiskDHI(1.0f0*GalaxyProfiles.defaultunits.length,1.0f0*GalaxyProfiles.defaultunits.mass) isa ExponentialDisk{Float32}
        @test ExponentialDiskDHI(1.0f0*GalaxyProfiles.defaultunits.length,1.0f0*GalaxyProfiles.defaultunits.mass,1f6*GalaxyProfiles.defaultunits.surfacedensity) isa ExponentialDisk{Float32}

    end
end
