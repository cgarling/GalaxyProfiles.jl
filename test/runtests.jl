using GalaxyProfiles
using Test

# first run regular tests, without units
tests = ["general_isothermal.jl","exponential_disk.jl","nfw.jl"]
for i in tests
    include(i)
end

# now run tests for Unitful API
import Unitful as u
import UnitfulAstro as ua
tests_units = ["general_isothermal_units.jl","exponential_disk_units.jl"]

for i in tests_units
    include(i)
end
