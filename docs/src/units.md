# [Integration with Unitful.jl](@id units)
`GalaxyProfiles.jl` provides integration with [`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl) and [`UnitfulAstro.jl`](https://github.com/JuliaAstro/UnitfulAstro.jl) through [`Requires.jl`](https://github.com/JuliaPackaging/Requires.jl). When both `Unitful` and `UnitfulAstro` are imported, `GalaxyProfiles.jl/src/units.jl` will be included, which will overload many of the default [methods](@ref methods) and [constructors](@ref types) to work with `Unitful` quantities.

## Type Constructors
Constructors for our defined types will accept Unitful quantities with correct dimensionality. So as to not explicitly depend on `Unitful` and `UnitfulAstro`, units are never stored in types. Instead, they are converted to [default units](@ref defaultunits) and stored in the types as `Real`s.

```@repl units
import Unitful as u
import UnitfulAstro as ua
using GalaxyProfiles

@info "Correct dimensions, converted and stored internally as Float64." # hide
d = GeneralIsothermal(1.0 * ua.Msun / ua.kpc^3, 1.0 * ua.kpc, 2.5)

@info "Incorrect dimensions on first argument, will error." # hide
d = GeneralIsothermal(1.0 * ua.Msun / ua.kpc^2, 1.0 * ua.kpc, 2.5)
```

## [Default Units](@id defaultunits)

The default units used by type constructors and methods are defined in the submodule `GalaxyProfiles.defaultunits` (found at the top of `src/units.jl`) and are generally named after the quantity that they represent. These include
 - mass: `UnitfulAstro.Msun`
 - ∇mass: `UnitfulAstro.Msun / UnitfulAstro.kpc`
 - density: `UnitfulAstro.Msun / UnitfulAstro.kpc^3`
 - surfacedensity: `UnitfulAstro.Msun / UnitfulAstro.kpc^2`
 - length: `UnitfulAstro.kpc`
 - velocity: `Unitful.km / Unitful.s`
 - time: `Unitful.yr`
 - Φunit: `Unitful.km^2 / Unitful.s^2`
 - ∇Φunit: `Unitful.km^2 / Unitful.s^2 / UnitfulAstro.kpc`
 - ∇∇Φunit: `Unitful.km^2 / Unitful.s^2 / UnitfulAstro.kpc^2`

Accompanying dimensions for proper dispatch are also defined in `src/units.jl`, including `SurfaceDensity, ∇ρdimension, Φdimension, ∇∇Φdimension, ∇mdimension`. These function like so:

```@repl units
@info "Get Unitful package extension, if on version of Julia that supports it" # hide
if isdefined(Base, :get_extension)
    ext = Base.get_extension(GalaxyProfiles, :GalaxyProfilesUnitfulExt)
else
    ext = GalaxyProfiles.GalaxyProfilesUnitfulExt
end
1 * ua.Msun / ua.pc^2 isa ext.SurfaceDensity
ua.Msun / ua.pc^2 isa ext.SurfaceDensityUnits
```

See also the documentation for [`Unitful.@derived_dimension`](https://painterqubits.github.io/Unitful.jl/stable/newunits/#Unitful.@derived_dimension).

## Methods

Methods defined on `Real` inputs are not overloaded to return Unitful quantities. This is to ensure identical behaviors of these methods regardless of whether Unitful integration is active or not.

```@repl units
ρ(GeneralIsothermal(1.0 * ua.Msun / ua.kpc^3, 1.0 * ua.kpc, 2.5), 1.0)
```

**It is therefore safest when constructing types with Unitful inputs to also call methods with Unitful inputs.** Otherwise, you can accidentally end up mixing units in ways that are hard to reason about. See the section on [units warning](@ref unitswarning).

Methods that only take a [`GalaxyProfiles.AbstractMassProfile`](@ref) instance as input generally have one additional method allowing for a conversion of that quantity. For example,

```@repl units
d = ExponentialDisk(1 * ua.Msun / ua.pc^2, 100 * ua.pc)
Mtot(d) # default method
Mtot(ua.Msun, d)
Mtot(u.kg, d)
```

Methods that take a [`GalaxyProfiles.AbstractMassProfile`](@ref) and a `Real` will have three additional methods allowing for Unitful integrations. For example,

```@repl units
d = GeneralIsothermal(1.0 * ua.Msun / ua.kpc^3, 1.0 * ua.kpc, 2.5)
ρ(d, 1.0) # default method
ρ(ua.Msun/ua.pc^3, d, 1.0) # for real argument and a unit conversion of result
ρ(d, 1.0*ua.kpc) # for a Unitful argument
ρ(ua.Msun/ua.pc^3, d, 1.0*ua.kpc) # for Unitful argument and result conversion
```

### [Units Warning](@id unitswarning)
For example, say that I define the scale radius in parsecs when constructing a type. This will be converted and stored in the type as a plain float type with an implicit unit of kpc (see above section on [default units](@ref defaultunits)). If I then call a method with a plain float radius argument that I expect to be in the same units as the scale radius I provided, the result will be in the wrong units. **After constructing a type with Unitful quantities, all methods called on that type with plain float inputs must be in the [default units](@ref defaultunits) that are expected**. When calling methods with basic numeric inputs, you are responsible for managing the units properly. A solution is to provide units on inputs when calling methods to ensure all units are properly converted. Calls like this will return Unitful quantities. 

```@repl units
@info "Unitful arguments are converted to basic numeric types."
d = GeneralIsothermal(1.0 * ua.Msun / ua.kpc^3, 1.0 * ua.pc, 2.5)
rs = 1.0 # I might expect this to be in parsecs, like the scale radius I provided ...
ρ(d, rs) # But this assumes that rs is in the default length unit, which is kpc.
rs = 1.0 * ua.pc # The safer way to do this is provide rs with units.
ρ(d, rs)
```