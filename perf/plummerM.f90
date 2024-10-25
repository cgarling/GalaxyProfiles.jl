! Mass enclosed within radius r for a Plummer profile with total mass M and scale radius a.
! In Julia this is
! a * M * r^3 * sqrt(1 + (r/a)^2) / (a^2 + r^2)^2
! When using REAL (single-precision) variables, the runtimes are as follows: 
! Compiling as gfortran plummerM.f90 -o plummerM  gives runtime of ~340 microseconds.
! Compiling as gfortran plummerM.f90 -o plummerM  -O2 gives runtime of ~260 microseconds. This seems to be where a lot of guides recommend setting the optimization level.
! Compiling as gfortran plummerM.f90 -o plummerM -O3 gives runtime of ~90 microseconds.
! Compiling as gfortran plummerM.f90 -o plummerM -Ofast gives runtime of ~80 microseconds. Adding -funroll-loops doesn't change it. Using -Ofast gives similar timings but actualy a different result, presumably due to associativity of the loop sum.
! See commented section at end of file for timings from Julia implementation. 
! AH, part of it is that fortran REAL is single-precision by default. Need to declar DOUBLE PRECISION to be equivalent to Julia Float64; Julia's Float32 and Fortran REAL give the same result.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! When running in double precision, the runtimes are as follows:
! Compiling as gfortran plummerM.f90 -o plummerM  gives runtime of ~430 microseconds, ~40% slower than Julia. 
! Compiling as gfortran plummerM.f90 -o plummerM -O2 gives runtime of ~400 microseconds, ~33% slower than Julia. -O3 is the same actually.
! Compiling as gfortran plummerM.f90 -o plummerM -Ofast gives runtime of ~140 microseconds, about half the time of the default Julia implementation, but slower than LoopVectorization. Adding -funroll-loops doesn't change it. Using -Ofast gives similar timings but actualy a different result, presumably due to associativity of the loop sum.

PROGRAM  plummerM
   IMPLICIT  NONE

   ! REAL     :: r, M, a, ans ! evaluation radius, mass, scale radius, the result, and the execution time.
   DOUBLE PRECISION     :: r, M, a, ans ! evaluation radius, mass, scale radius, the result
   REAL     :: startT, endT, execTime ! for keeping track of time
   INTEGER  :: i ! Loop indexing.
   WRITE(*,*) "Enter evaluation radius, total mass, and Plummer scale radius:"
   ! READ(*,*)  r, M, a
   M = 100.0
   a = 20.0
   ans = 0.0

   call cpu_time(startT) ! Start timer
   ! This is too short to profile really. Make a loop and do a simple sum. 
   ! ans = a * M * r**3 * SQRT(1 + (r/a)**2) / (a**2 + r**2)**2
   DO i = 1, 10**5, 1
      r = i * 1E-3
      ans = ans + a * M * r**3 * SQRT(1 + (r/a)**2) / (a**2 + r**2)**2
   END DO
   call cpu_time(endT)
   execTime = endT - startT
   WRITE(*,*) "Enclosed Mass = ", ans
   WRITE(*,*) 'Elapsed time in seconds: ', execTime

END PROGRAM  plummerM

! Equivalent Julia code
! In single precision,
! rr2 = (1:10^5) .* 1f-3 |> collect
! sum( (M(Plummer(100.0f0,20.0f0),i) for i in rr2))
! @benchmark sum( (M(Plummer($100.0f0,$20.0f0),i) for i in $rr2))
! takes ~140 microseconds. 
! In double precision,
! rr = (1:10^5) .* 1e-3 |> collect
! sum( (M(Plummer(100.0,20.0),i) for i in rr))
! @benchmark sum( (M(Plummer($100.0,$20.0),i) for i in $rr))
! runs in ~300 microseconds
! import LoopVectorization: @turbo
! function psum(arr, a, M)
!   result = zero(typeof(a))
!   # This does not give any speedup, but @turbo does
!   # @inbounds @simd for i in eachindex(arr)
!   # for i in eachindex(arr)
!   @turbo for i in eachindex(arr)
!     r = arr[i]
!     result += a * M * r^3 * sqrt(1 + (r/a)^2) / (a^2 + r^2)^2
!   end
!   return result
! end
! This will also work;
! function psum2(d::Plummer, rr::AbstractArray{<:Real})
!     result = zero(eltype(rr))
!     @turbo for i in eachindex(rr)
!         result += M(d, rr[i])
!     end
!     return result
! end
! @benchmark psum($rr, $20.0, $100.0)
! Runs in about 300 microseconds in double precision with @inbounds @simd, which is the same as an un-annotated loop. 
! Runs in about 62 microseconds in single precision with @inbounds @simd, and 192 microseconds with a simple for loop. 
! Runs in about 100 microseconds in double precision with LoopVectorization, and 28 microseconds in single precision. 
