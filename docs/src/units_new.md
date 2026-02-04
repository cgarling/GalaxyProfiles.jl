# [Integration with Unitful.jl](@id units)

`GalaxyProfiles.jl` provides seamless integration with [`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl) and [`UnitfulAstro.jl`](https://github.com/JuliaAstro/UnitfulAstro.jl). Mass profile types can store `Unitful.Quantity` values directly in their fields, and all methods work transparently with these quantities, automatically propagating units through calculations.

## Type Constructors

Constructors for mass profile types accept both plain `Real` numbers and `Unitful` quantities with correct dimensionality. **Units are preserved** and stored directly in the struct fields, ensuring dimensional correctness throughout all calculations.

```@repl units
import Unitful as u
import UnitfulAstro as ua
using GalaxyProfiles

@info "Creating Plummer with Unitful quantities - units are preserved" # hide
d = Plummer(1e8 * ua.Msun, 0.5 * ua.kpc)
d.M  # Returns 1.0e8 M⊙ 
d.a  # Returns 0.5 kpc

@info "Creating Plummer with plain numbers" # hide
d_plain = Plummer(1.0, 1.0)  # Plain Float64 values
d_plain.M  # Returns 1.0 (dimensionless)
```

**Key Difference from Previous Versions:**
- **New behavior:** Units are preserved in struct fields and propagate through all calculations
- **Old behavior:** Units were stripped during construction and stored as plain floats

## Working with Unitful Quantities

When you create a mass profile with Unitful quantities, all methods automatically return results with appropriate units:

```@repl units
d = Plummer(1e8 * ua.Msun, 0.5 * ua.kpc)

# Density at a radius
ρ(d, 0.3 * ua.kpc)  # Returns density with units M⊙ kpc^-3

# Surface density  
Σ(d, 0.3 * ua.kpc)  # Returns surface density with units M⊙ kpc^-2

# Enclosed mass
M(d, 1.0 * ua.kpc)  # Returns mass with units M⊙

# Sampling functions return radii with units
r_samples = sample2D_r(d, 10)  # Returns array of radii with units kpc
```

## Mixed Units

You can mix different (but compatible) units - Julia's Unitful will handle conversions automatically:

```@repl units
d1 = Plummer(1e8 * ua.Msun, 0.5 * ua.kpc)
d2 = Plummer(1e8 * ua.Msun, 500.0 * ua.pc)  # Using parsecs instead

# These give equivalent results at the same physical radius
ρ(d1, 0.5 * ua.kpc)
ρ(d2, 500.0 * ua.pc)
```

## Backwards Compatibility

Plain numeric types (without units) continue to work exactly as before:

```@repl units
# Plain numbers - no units
d = Plummer(1.0, 1.0)
ρ(d, 1.0)  # Returns plain Float64

# When using plain numbers, default units are assumed:
# [M] = M⊙, [r, a] = kpc
# This is important for gravitational quantities like Vcirc, Vesc, Φ
```

## Default Units for Plain Numbers

When working with plain `Real` numbers (no Unitful quantities), the following default units are assumed:

 - mass: `M⊙` (solar masses)
 - length/radius: `kpc` (kiloparsecs)
 - density: `M⊙ kpc^-3`
 - surface density: `M⊙ kpc^-2`
 - velocity: `km s^-1`
 - time: `yr` (years)

These defaults are important for physical quantities that involve the gravitational constant `G`, such as [`Vcirc`](@ref), [`Vesc`](@ref), [`Φ`](@ref), [`∇Φ`](@ref), and [`∇∇Φ`](@ref).

## Dimensional Safety

Using Unitful quantities provides automatic dimensional analysis - operations with incompatible dimensions will produce clear error messages:

```@repl units
d = Plummer(1e8 * ua.Msun, 0.5 * ua.kpc)

# This will error - can't use a mass where a length is expected:
# ρ(d, 1.0 * ua.Msun)  # DimensionError!

# This will error - mixing incompatible dimensions:
# d2 = Plummer(1e8 * ua.kpc, 0.5 * ua.Msun)  # Wrong parameter types!
```

## Type Parameters

Mass profile types use parametric types to support both plain numbers and Unitful quantities:

```julia
# Type signatures allow flexibility:
struct Plummer{TM, Ta, T} <: AbstractDensity{T}
    M::TM  # Can be Float64 or Unitful.Quantity
    a::Ta  # Can be Float64 or Unitful.Quantity
end

# Examples:
Plummer{Float64, Float64, Float64}(1.0, 1.0)
Plummer{Quantity{Float64, 𝐌, ...}, Quantity{Float64, 𝐋, ...}, ...}(1e8*ua.Msun, 0.5*ua.kpc)
```

The type parameter `T` represents the result type of operations involving the fields, computed as the product of the field types.
