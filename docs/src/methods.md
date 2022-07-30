# [Defined Methods](@id methods)
These common methods are defined for most density profiles, when possible.

```@docs
ρ(::GalaxyProfiles.AbstractDensity,::Real)
invρ(::GalaxyProfiles.AbstractDensity,::Real,::Tuple)
∇ρ(::GalaxyProfiles.AbstractDensity,::Real)
ρmean(::GalaxyProfiles.AbstractDensity,::Real)
invρmean(::GalaxyProfiles.AbstractDensity,::Real,::Tuple)
Σ(::GalaxyProfiles.AbstractMassProfile,::Real)
invΣ(::GalaxyProfiles.AbstractMassProfile,::Real)
∇Σ(::GalaxyProfiles.AbstractMassProfile,::Real)
Σmean(::GalaxyProfiles.AbstractMassProfile,::Real)
M(::GalaxyProfiles.AbstractDensity,::Real)
invM(::GalaxyProfiles.AbstractDensity,::Real)
∇M(::GalaxyProfiles.AbstractDensity,::Real)
Mtot(::GalaxyProfiles.AbstractDensity)
Mproj(::GalaxyProfiles.AbstractMassProfile,::Real)
∇Mproj(::GalaxyProfiles.AbstractMassProfile,::Real)
invMproj(::GalaxyProfiles.AbstractMassProfile,::Real)
cdf
ccdf
quantile
cquantile
Vcirc(::GalaxyProfiles.AbstractDensity,::Real)
Vesc(::GalaxyProfiles.AbstractDensity,::Real)
Φ(::GalaxyProfiles.AbstractDensity,::Real)
∇Φ(::GalaxyProfiles.AbstractDensity,::Real)
∇∇Φ(::GalaxyProfiles.AbstractDensity,::Real)
```