"""
    ρ(d::AbstractDensity,r::Real)
    ρ(d::AbstractDensity,r::Unitful.Quantity)

Evaluate the density of the `AbstractDensity` profile at radius `r`.
"""
ρ()

"""
    Σ(d::AbstractDensity,r::Real)
    Σ(d::AbstractDensity,r::Unitful.Quantity)
    Σ(d::AbstractSurfaceDensity,r::Real)
    Σ(d::AbstractSurfaceDensity,r::Unitful.Quantity)

Evaluate the surface density of the `AbstractDensity` or `AbstractSurfaceDensity` profile at radius `r`.
"""
Σ()
