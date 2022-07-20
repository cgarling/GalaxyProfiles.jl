################################################################
### Default implementations for random sampling

"""
    rand!([rng::Random.AbstractRNG=Random.default_rng()], d::AbstractMassProfile, x::AbstractArray{<:Real})

Draw multiple samples from `d` and write them to `x` in place. By default, falls back to equivalent of `quantile(d,rand(rng))`. 
"""
function rand!(rng::AbstractRNG,d::AbstractMassProfile,x::AbstractArray{<:Real})
    rand!(rng,x)
    x .= quantile.(d,x)
end
rand!(d::AbstractMassProfile,x::AbstractArray{<:Real}) = rand!(default_rng(),d,x)
"""
    rand([rng=Random.default_rng()], d::AbstractMassProfile, [dims...])

Sample random points from `d` with shape `[dims...]`. By default, falls back to equivalent of `quantile(d,rand(rng))`.
"""
rand(rng::AbstractRNG,d::AbstractMassProfile,dims::Int...) = quantile.(d,rand(rng,dims...))
rand(d::AbstractMassProfile,dims::Int...) = rand(default_rng(),d,dims...)
rand(rng::AbstractRNG,d::AbstractMassProfile,dims::Dims) = quantile.(d,rand(rng,dims...))
rand(d::AbstractMassProfile,dims::Dims) = rand(default_rng(),d,dims)
rand(rng::AbstractRNG,d::AbstractMassProfile) = quantile(d,rand(rng))
rand(d::AbstractMassProfile) = rand(default_rng(),d)
