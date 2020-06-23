"""
$(TYPEDSIGNATURES)

```jldoctest
julia> random_design!(zeros(10, 3), (Uniform(2, 3), DiscreteUniform(-1, 5), Uniform(5, 10)), 10)
10×3 Array{Float64,2}:
 2.04922  -1.0  6.44508
 2.59117   3.0  5.57183
 2.77148   5.0  8.72759
 2.25659   1.0  9.75865
 2.64968   4.0  5.72029
 2.31523   5.0  8.48192
 2.86526   5.0  5.42517
 2.44753   4.0  6.75815
 2.08284   3.0  8.94993
 2.81201  -1.0  5.89413

```
"""
function random_design!(design::Array{Float64, 2}, distributions::Tuple, n::Int)
    for d in 1:size(distributions, 1)
        design[:, d] .= rand(distributions[d], n)
    end

    return design
end
