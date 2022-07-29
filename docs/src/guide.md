# [Defining a New Profile](@id guide)
If your profile describes a 3D density distribution, create a new file under the `src/densities` directory. If only 2D quantities will be defined, create a new file under `src/surface_densities` instead. Let's say I want to define a special type for the singular isothermal sphere, rather than [`SIS`](@ref) that creates a [`GeneralIsothermal`](@ref) instance. I could create `src/sis.jl` and add `include("sis.jl")` to `src/densities.jl` to make sure the file is loaded. Then, in `src/sis.jl` I could define a singular isothermal sphere as

```@setup guide
import GalaxyProfiles:AbstractDensity
```

```@example guide
"""
    SingularIsothermalSphere(ρ0::Real,rs::Real)

Singular isothermal sphere profile. Fields are `ρ0,rs`. 
"""
struct SingularIsothermalSphere{T<:Real} <: AbstractDensity
    ρ0::T
    rs::T
end
nothing # hide
```

You should add some info about the profile to the docstring. If you are creating a 2D surface density distribution, you should subtype `AbstractSurfaceDensity`. You should then define additional constructors for your type. To allow the call signature we listed in the docstring above, define

```@example guide
SingularIsothermalSphere(ρ0::Real,rs::Real) = SingularIsothermalSphere(promote(ρ0,rs)...)
nothing # hide
```

such that we can do

```@example guide
SingularIsothermalSphere(1.0,1)
```

It is also common to define constructors that will take the total mass and scale radius as inputs and compute the characteristic density; in this example,`ρ0`. See the implementation of [`GeneralIsothermal`](@ref) for an example of how to do this.

Now you should define as many methods from the [Defined Methods](@ref methods) section as you need, starting with methods to retrieve the parameters,

```@example guide
params(d::SingularIsothermalSphere) = (d.ρ0,d.rs)
scale_radius(d::SingularIsothermalSphere) = d.rs
nothing # hide
```

and then evaluation methods like

```@example guide
function ρ(d::SingularIsothermalSphere,r::Real)
    ρ0, rs = params(d)
    ρ0 * (r/rs)^-2
end
nothing # hide
```

We avoid accessing the fields of the types directly and use `params` and `scale_radius` instead so that internal fields of the types can be refactored later if necessary without having to redefine all the accompanying methods.

When you are done writing your methods, you should define methods for your type that allow for `Unitful` quantities in the `src/units.jl` file. You should be able to follow the examples there with little problem. You should write tests under `GalaxyProfiles.jl/test` to validate the behavior of your type and methods. It is recommended that you make a new file, e.g. `test/mytype.jl`, then add `include("mytype.jl")` to the `test/runtests.jl` file. You should then edit `src/GalaxyProfiles.jl` to export your new type, in this example `SingularIsothermalSphere`. 