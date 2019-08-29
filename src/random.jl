"""
$(TYPEDSIGNATURES)

```jldoctest
julia> random_design((Uniform(2, 3), DiscreteUniform(-1, 5), Uniform(5, 10)), 10)
10×3 Array{Real,2}:
 2.04922   5  8.85741
 2.25659  -1  6.57617
 2.86526   4  5.41422
 2.81201  -1  5.40944
 2.92412   1  6.22902
 2.20051   5  8.80749
 2.69459   3  6.34626
 2.28902   2  8.72759
 2.95173   0  5.42517
 2.35163   4  5.89413

```
"""
function random_design(distributions::Tuple, size::Int)
    permutedims(foldl((x, y) -> hcat(collect(x),
                                     collect(y)),
                      (rand(d) for d in distributions) for i = 1:size))
end
