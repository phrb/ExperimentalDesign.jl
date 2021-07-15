"""
$(TYPEDSIGNATURES)

```jldoctest
julia> random_design!((Uniform(2, 3), DiscreteUniform(-1, 5), Uniform(5, 10)), 10)
10×3 Matrix{Float64}:
 2.04922  -1.0  8.21161
 2.59117   3.0  5.40944
 2.77148  -1.0  9.62061
 2.25659   4.0  5.79254
 2.64968   2.0  6.22902
 2.31523   5.0  6.00256
 2.86526   2.0  9.21818
 2.44753   0.0  8.80749
 2.08284  -1.0  8.47297
 2.81201  -1.0  6.11722

```
"""
function random_design!(distributions::Tuple, n::Int)
    return hcat([rand(distributions[d], n) for d in 1:size(distributions, 1)]...)
end
