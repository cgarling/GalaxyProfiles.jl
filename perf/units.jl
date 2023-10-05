# using GalaxyProfiles
import Unitful as u
import UnitfulAstro as ua
import BenchmarkTools:@btime,@benchmark

struct ExponentialDisk{T<:Real}
    Σ0::T
    rs::T
end
Base.Broadcast.broadcastable(m::ExponentialDisk) = Ref(m)
Σ(d::ExponentialDisk,r::Real) = d.Σ0 * exp(-r/d.rs)
function Mproj(d::ExponentialDisk,r::Real)
    isinf(r) && (return 2π * d.Σ0 * d.rs^2)
    ee = exp(-r/d.rs)
    2π * d.Σ0 * d.rs^2 * (1 - ee*(r/d.rs + 1))
end
struct ExponentialDiskUnits{T<:Number,S<:Number}
    Σ0::T
    rs::S
end
Base.Broadcast.broadcastable(m::ExponentialDiskUnits) = Ref(m)
Σ(d::ExponentialDiskUnits,r::Number) = d.Σ0 * exp(-r/d.rs)
function Mproj(d::ExponentialDiskUnits,r::Number)
    isinf(r) && (return 2π * d.Σ0 * d.rs^2)
    ee = exp(-r/d.rs)
    2π * d.Σ0 * d.rs^2 * (1 - ee*(r/d.rs + 1))
end
# this actually works fine; might be worth implementing
# struct ExponentialDiskUnits{T<:Real,S<:u.Quantity{T,u.𝐌/u.𝐋^2},V<:u.Quantity{T,u.𝐋}}
#     Σ0::S
#     rs::V
# end
# the T is not technically necessary but nice to have I guess ? ^^^
# See discussion
# https://discourse.julialang.org/t/how-to-properly-use-unitful/40295/3
Σ(d::ExponentialDiskUnits,r::u.Length)= d.Σ0 * exp(-r/d.rs)

Base.Broadcast.broadcastable(m::ExponentialDiskUnits) = Ref(m)
function testbasic()
    Σ0_nou = 1e6
    Σ0_u = 1.0 * ua.Msun / ua.pc^2 |> ua.Msun/ua.kpc^2
    rs_nou = 10.0
    rs_u = 10.0 * ua.kpc
    r_nou = 1.0
    r_u = 1.0 * ua.kpc
    r_u2 = 1.0e3  * ua.pc
    println("Struct Creation; first no units, then units")
    @btime ExponentialDisk($Σ0_nou,$rs_nou)
    @btime ExponentialDiskUnits($Σ0_u,$rs_u)
    d_nou = ExponentialDisk(Σ0_nou,rs_nou)
    d_u = ExponentialDiskUnits(Σ0_u,rs_u)
    println("Surface Density Evaluation; first no units, then units")
    @btime Σ($d_nou,$r_nou)
    @btime Σ($d_u,$r_u)
    println("Projected Mass Evaluation; first no units, then units")
    @btime Mproj($d_nou,$r_nou)
    @btime Mproj($d_u,$r_u)
    println("Check semantics of units when providing mismatching arguments")
    println(Σ(d_u,r_u))
    println(Σ(d_u,r_u2))
    println("Look at performance of broadcasting")
    rvec_nou = rand(1000)
    rvec_u = copy(rvec_nou)*ua.kpc
    @btime Σ.($d_nou,$rvec_nou)
    @btime Σ.($d_u,$rvec_u)

    nothing
end
testbasic()
# looks like there's no real drawback to using units. Should maybe implement them? 
