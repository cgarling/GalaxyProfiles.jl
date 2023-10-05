# Looking at vectorization and performance of Float64 vs Float32
import BenchmarkTools: @benchmark
using GalaxyProfiles
import LoopVectorization: @turbo

mass = 100.0
a  = 20.0
rsingle = (1:10^5) .* 1f-3
rdouble = (1:10^5) .* 1e-3
benchmark_time = 0.5 # seconds

# What's the takeaway?
# The M(::Plummer, ::Real) interface has no overhead.
# In single precision (Float32) we get SIMD whereas we don't in double precision. This makes a significant speed different 300 μs in double |> 60 in single.
# LoopVectorization does not like the M(::Plummer, ::Real) interface and will break. In single precision LV gives a 2x speedup and in double precision it gives a 3x speedup. All of these are competitive / better than the equivalent Fortran timings for the program in plummerM.f90.
# Overall looking very good but making me question whether I should work in single precision sometimes ... 

# Msum computes sum( r-> M(Plummer(mass,a), r), arr ) using different methods
# This method does not use the M(Plummer(mass,a), r) interface, but rather the equation directly
function Msum_base(arr, a, mass)
    result = zero(typeof(a))
    for i in eachindex(arr)
        r = arr[i]
        result += a * mass * r^3 * sqrt(1 + (r/a)^2) / (a^2 + r^2)^2
    end
    return result
end
function Msum_simd(arr, a, mass)
    result = zero(typeof(a))
    @inbounds @simd for i in eachindex(arr)
        r = arr[i]
        result += a * mass * r^3 * sqrt(1 + (r/a)^2) / (a^2 + r^2)^2
    end
    return result
end
function Msum_turbo(arr, a, mass)
    result = zero(typeof(a))
    @turbo for i in eachindex(arr)
        r = arr[i]
        result += a * mass * r^3 * sqrt(1 + (r/a)^2) / (a^2 + r^2)^2
    end
    return result
end

println("Enclosed Mass of Plummer Profile: Full Equation")
println("=================================================================")
println("Julia `sum`, full equation, single precision: ")
# sum( r-> Float32(a) * Float32(mass) * r^3 * sqrt(1 + (r/Float32(a))^2) / (Float32(a)^2 + r^2)^2, rsingle )
# b = @benchmark sum( r-> Float32($a) * Float32($mass) * r^3 * sqrt(1 + (r/Float32($a))^2) / (Float32($a)^2 + r^2)^2, $rsingle ) seconds=0.5
let a = Float32(a), mass = Float32(mass)
    b = @benchmark sum( r-> $a * $mass * r^3 * sqrt(1 + (r/$a)^2) / ($a^2 + r^2)^2, $rsingle ) seconds=0.5
    show(stdout, "text/plain", b)
end
println("")
println("=================================================================")  
println("Basic loop, full equation, single precision: ")
b = @benchmark Msum_base( $rsingle, $Float32(a), $Float32(mass) ) seconds=benchmark_time
show(stdout, "text/plain", b)
println("")
println("=================================================================")         
println("@inbounds @simd loop, full equation, single precision: ")
b = @benchmark Msum_simd( $rsingle, $Float32(a), $Float32(mass) ) seconds=benchmark_time
show(stdout, "text/plain", b)
println("")
println("=================================================================")
println("LoopVectorization @turbo loop, full equation, single precision: ")
b = @benchmark Msum_turbo( $rsingle, $Float32(a), $Float32(mass) ) seconds=benchmark_time
show(stdout, "text/plain", b)
println("")
println("=================================================================")

println("Julia `sum`, full equation, double precision: ")
# sum( r-> a * mass * r^3 * sqrt(1 + (r/a)^2) / (a^2 + r^2)^2, rdouble )
b = @benchmark sum( r-> $a * $mass * r^3 * sqrt(1 + (r/$a)^2) / ($a^2 + r^2)^2, $rsingle ) seconds=0.5
show(stdout, "text/plain", b)
println("")
println("=================================================================")
println("Basic loop, full equation, double precision: ")
b = @benchmark Msum_base( $rdouble, $a, $mass ) seconds=benchmark_time
show(stdout, "text/plain", b)
println("")
println("=================================================================")
println("@inbounds @simd loop, full equation, double precision: ")
b = @benchmark Msum_simd( $rdouble, $a, $mass ) seconds=benchmark_time
show(stdout, "text/plain", b)
println("")
println("=================================================================")
println("LoopVectorization @turbo loop, full equation, double precision: ")
b = @benchmark Msum_turbo( $rdouble, $a, $mass ) seconds=benchmark_time
show(stdout, "text/plain", b)
println("")
println("=================================================================")

# Now these use the M(Plummer(mass,a), r) interface
function Msum2_base(arr, a, mass)
    result = zero(typeof(a))
    d = Plummer(mass, a)
    for i in eachindex(arr)
        r = arr[i]
        result += M(d, r)
    end
    return result
end
function Msum2_simd(arr, a, mass)
    result = zero(typeof(a))
    d = Plummer(mass, a)
    @inbounds @simd for i in eachindex(arr)
        r = arr[i]
        result += M(d, r)
    end
    return result
end
# LoopVectorization does not like the custom struct `d` but since we aren't
# indexing into it or looping over it I don't understand what the problem is really. 
function Msum2_turbo(arr, a, mass)
    result = zero(typeof(a))
    d = Plummer(mass, a)
    @turbo for i in eachindex(arr)
        r = arr[i]
        result += M(d, r)
    end
    return result
end

println("Enclosed Mass of Plummer Profile: Plummer Type Interface")
println("=================================================================")
println("Julia `sum`, Plummer interface, single precision: ")
# sum( r-> M(Plummer(Float32($mass),Float32($a)), r), rsingle )
# b = @benchmark sum( r-> M(Plummer(Float32($mass),Float32($a)), r), $rsingle ) seconds=0.5
let d = Plummer(Float32(mass), Float32(a))
    b = @benchmark sum( r-> M($d, r), $rsingle ) seconds=0.5
    show(stdout, "text/plain", b)
end
println("")
println("=================================================================")
println("Basic loop, Plummer interface, single precision: ")
b = @benchmark Msum2_base( $rsingle, $Float32(a), $Float32(mass) ) seconds=benchmark_time
show(stdout, "text/plain", b)
println("")
println("=================================================================")         
println("@inbounds @simd loop, Plummer interface, single precision: ")
b = @benchmark Msum2_simd( $rsingle, $Float32(a), $Float32(mass) ) seconds=benchmark_time
show(stdout, "text/plain", b)
println("")
println("=================================================================")
println("LoopVectorization @turbo loop, Plummer interface, single precision: ")
b = @benchmark Msum2_turbo( $rsingle, $Float32(a), $Float32(mass) ) seconds=benchmark_time
show(stdout, "text/plain", b)
println("")
println("=================================================================")
println("Julia `sum`, Plummer interface, double precision: ")
# sum( r-> M(Plummer(mass,a), r), rdouble )
# b = @benchmark sum( r-> M(Plummer($mass,$a), r), $rdouble ) seconds=0.5
let d = Plummer(mass, a)
    b = @benchmark sum( r-> M($d, r), $rdouble ) seconds=0.5
    show(stdout, "text/plain", b)
end
println("")
println("=================================================================")
println("Basic loop, Plummer interface, double precision: ")
b = @benchmark Msum2_base( $rdouble, $a, $mass ) seconds=benchmark_time
show(stdout, "text/plain", b)
println("")
println("=================================================================")
println("@inbounds @simd loop, Plummer interface, double precision: ")
b = @benchmark Msum2_simd( $rdouble, $a, $mass ) seconds=benchmark_time
show(stdout, "text/plain", b)
println("")
println("=================================================================")
println("LoopVectorization @turbo loop, Plummer interface, double precision: ")
b = @benchmark Msum2_turbo( $rdouble, $a, $mass ) seconds=benchmark_time
show(stdout, "text/plain", b)
println("")
println("=================================================================")




