# Defined Types
## Abstract Types
The highest level type defined in this package is the abstract type [`AbstractMassProfile`](@ref GalaxyProfiles.AbstractMassProfile), which all other profiles are subtyped from. Below these are [`AbstractDensity`](@ref GalaxyProfiles.AbstractDensity) and [`AbstractSurfaceDensity`](@ref GalaxyProfiles.AbstractSurfaceDensity); the latter is defined only as a surface density and will have no 3D quantities defined (e.g., [`ρ`](@ref)), while the former represents 3D density profiles.
```@docs
GalaxyProfiles.AbstractMassProfile
GalaxyProfiles.AbstractDensity
GalaxyProfiles.AbstractSurfaceDensity
```

## Concrete Types
The following concrete types, representing specific density profiles, are currently implemented:

```@docs
ExponentialDisk
GeneralIsothermal
```

## Retrieving Parameters
The parameters that define these types can be retrieved with [`params`](@ref) and [`scale_radius`](@ref)
```@docs
GalaxyProfiles.params(::GalaxyProfiles.AbstractMassProfile)
scale_radius(::GalaxyProfiles.AbstractMassProfile)
```

## Convenience Constructors
We also provide some convenience constructors for other types such as the [`singular isothermal sphere`](@ref SIS), which returns an instance of [`GeneralIsothermal`](@ref) with `α=2`.
```@docs
SIS
```