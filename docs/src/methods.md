# [Defined Methods](@id methods)

## Common Methods
These common methods are defined for most density profiles, when possible.

```@docs
ρ(::GalaxyProfiles.AbstractDensity,::Real)
invρ(::GalaxyProfiles.AbstractDensity,::T,::NTuple{2,S}) where {T<:Real,S<:Real}
∇ρ(::GalaxyProfiles.AbstractDensity,::Real)
ρmean(::GalaxyProfiles.AbstractDensity,::Real)
invρmean(::GalaxyProfiles.AbstractDensity,::T,::NTuple{2,S}) where {T<:Real,S<:Real}
Σ(::GalaxyProfiles.AbstractMassProfile,::Real)
invΣ(::GalaxyProfiles.AbstractMassProfile,::T,::NTuple{2,S}) where {T<:Real,S<:Real}
∇Σ(::GalaxyProfiles.AbstractMassProfile,::Real)
Σmean(::GalaxyProfiles.AbstractMassProfile,::Real)
M(::GalaxyProfiles.AbstractDensity,::Real)
invM(::GalaxyProfiles.AbstractDensity,::T,::NTuple{2,S}) where {T<:Real,S<:Real}
∇M(::GalaxyProfiles.AbstractDensity,::Real)
Mtot(::GalaxyProfiles.AbstractMassProfile)
Mproj(::GalaxyProfiles.AbstractMassProfile,::Real)
∇Mproj(::GalaxyProfiles.AbstractMassProfile,::Real)
invMproj(::GalaxyProfiles.AbstractMassProfile,::T,::NTuple{2,S}) where {T<:Real,S<:Real}
cdf2D
cdf3D
ccdf2D
ccdf3D
quantile2D
quantile3D
cquantile2D
cquantile3D(::GalaxyProfiles.AbstractDensity,::Real)
Vcirc(::GalaxyProfiles.AbstractDensity,::Real)
Vesc(::GalaxyProfiles.AbstractDensity,::Real)
Vmax(::GalaxyProfiles.AbstractDensity,::Real)
Φ(::GalaxyProfiles.AbstractDensity,::Real)
∇Φ(::GalaxyProfiles.AbstractDensity,::Real)
∇∇Φ(::GalaxyProfiles.AbstractDensity,::Real)
```

## Private Methods
The following methods are defined for convenience or internal use but are not exported.

```@docs
GalaxyProfiles.plummer_a_to_rh
GalaxyProfiles.plummer_angular_avalue
```

## Random Sampling
The following sampling methods for drawing positions from instantiated mass profiles are provided.

```@docs
sample2D_r!
sample2D_r
sample3D_r!
sample3D_r
```