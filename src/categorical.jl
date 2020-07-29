"""
$(TYPEDEF)

$(TYPEDFIELDS)

A simple wrapper for a `DiscreteUniform` distribution over non-numerical arrays.

"""
struct CategoricalFactor <: DiscreteUnivariateDistribution
    values::Array{Any, 1}
    distribution::DiscreteUniform
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> a = CategoricalFactor([:a, :b, 2, 1.0])
CategoricalFactor(
values: Any[:a, :b, 2, 1.0]
distribution: Distributions.DiscreteUniform(a=1, b=4)
)

```
"""
function CategoricalFactor(values::Array{Any, 1})
    CategoricalFactor(values, DiscreteUniform(1, length(values)))
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> rand(CategoricalFactor([:a, :b, 2, 1.0]), 6)
6-element Array{Any,1}:
  :b
 2
 1.0
 2
  :a
 2

```
"""
function rand(distribution::CategoricalFactor, n::Int = 1)
    [distribution.values[s] for s in rand(distribution.distribution, n)]
end
