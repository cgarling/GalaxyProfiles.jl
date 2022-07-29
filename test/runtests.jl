using GalaxyProfiles
using Test

tests = ["general_isothermal.jl"]

for i in tests
    include(i)
end
