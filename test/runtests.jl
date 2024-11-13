# using GalaxyProfiles
# using Test

# # first run regular tests, without units
# tests = ["general_isothermal.jl", "exponential_disk.jl", "nfw.jl", "corenfw_tests.jl", "plummer_tests.jl"]
# for i in tests
#     include(i)
# end

# # now run tests for Unitful API
# import Unitful as u
# import UnitfulAstro as ua
# tests_units = ["general_isothermal_units.jl", "exponential_disk_units.jl", "nfw_units.jl"]

# for i in tests_units
#     include(i)
# end

using Test, SafeTestsets

@testset "GalaxyProfiles.jl" verbose=true begin
    
    @safetestset "Common Fallback Methods" include("common_tests.jl")
    
    @testset "Densities" verbose=true begin
        @safetestset "NFW" include("densities/NFW_tests.jl")
        @safetestset "CoreNFW" include("densities/CoreNFW_tests.jl")
        @safetestset "GeneralIsothermal" include("densities/GeneralIsothermal_tests.jl")
        @safetestset "Plummer" include("densities/Plummer_tests.jl")
    end

    @testset "SurfaceDensities" verbose=true begin
        @safetestset "ExponentialDisk" include("surface_densities/ExponentialDisk_tests.jl")

    end

    @safetestset "Unitful.jl Extension" include("unitful_tests.jl")

end
