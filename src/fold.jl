"""
$(TYPEDSIGNATURES)

```jldoctest
julia> fold!(PlackettBurman(2))
8×3 DataFrame
 Row │ factor1  factor2  dummy1
     │ Int64    Int64    Int64
─────┼──────────────────────────
   1 │       1        1       1
   2 │       1       -1      -1
   3 │      -1       -1       1
   4 │      -1        1      -1
   5 │      -1       -1      -1
   6 │      -1        1       1
   7 │       1        1      -1
   8 │       1       -1       1

```
"""
function fold!(design::AbstractScreeningDesign)
    mirror_matrix = deepcopy(design.matrix)
    if mirror_matrix[1,1] isa Symbol
        for col in names(mirror_matrix)
            replace!(mirror_matrix[!, col], :low => :high, :high => :low)
        end
    else
        for col in names(mirror_matrix)
            replace!(mirror_matrix[!, col], 1 => -1, -1 => 1)
        end
    end
    append!(design.matrix, mirror_matrix)
end
