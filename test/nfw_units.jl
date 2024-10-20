# Just test constructors. Due to how we are overloading methods in units.jl via metaprogramming, if the overloaded methods work for one type they should work on others, as long as the call signature is the same. We are therefore not going to retest them all here, as they are mostly covered in general_isothermal_units.jl.
@testset "NFW Units" begin
    @testset "Float64" begin
        @test NFW(1.0*GalaxyProfiles.defaultunits.density, 1.0*GalaxyProfiles.defaultunits.length) isa NFW{Float64}
        @test NFW(1.0*GalaxyProfiles.defaultunits.density, 1.0*GalaxyProfiles.defaultunits.length) == NFW(1.0, 1.0)
    end
    @testset "Float32" begin
        @test NFW(1.0f0*GalaxyProfiles.defaultunits.density, 1*GalaxyProfiles.defaultunits.length) isa NFW{Float32}
        @test NFW(1.0f0*GalaxyProfiles.defaultunits.density, 1.0f0*GalaxyProfiles.defaultunits.length) == NFW(1.0f0, 1.0f0)
        @test NFW(1*GalaxyProfiles.defaultunits.density, 1.0f0*GalaxyProfiles.defaultunits.length) isa NFW{Float32}
    end
end
