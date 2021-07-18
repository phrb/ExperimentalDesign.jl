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
function fold!(design::PlackettBurman)
    fold_dataframe!(design.matrix)
end

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> fold!(FractionalFactorial2Level(@formula(y ~ a + b + a&b + c + a&c+ b&c + a&b&c)))
16×7 DataFrame
 Row │ a      b      c      a_b    a_c    b_c    a_b_c 
     │ Int64  Int64  Int64  Int64  Int64  Int64  Int64 
─────┼─────────────────────────────────────────────────
   1 │    -1     -1     -1      1      1      1     -1
   2 │     1     -1     -1     -1     -1      1      1
   3 │    -1      1     -1     -1      1     -1      1
   4 │     1      1     -1      1     -1     -1     -1
   5 │    -1     -1      1      1     -1     -1      1
   6 │     1     -1      1     -1      1     -1     -1
   7 │    -1      1      1     -1     -1      1     -1
   8 │     1      1      1      1      1      1      1
   9 │     1      1      1     -1     -1     -1      1
  10 │    -1      1      1      1      1     -1     -1
  11 │     1     -1      1      1     -1      1     -1
  12 │    -1     -1      1     -1      1      1      1
  13 │     1      1     -1     -1      1      1     -1
  14 │    -1      1     -1      1     -1      1      1
  15 │     1     -1     -1      1      1     -1      1
  16 │    -1     -1     -1     -1     -1     -1     -1

```
"""
function fold!(design::FractionalFactorial2Level)
    fold_dataframe!(design.matrix)
end

function fold_dataframe!(df::DataFrames.DataFrame)
    mirror_matrix = deepcopy(df)
    if mirror_matrix[1,1] isa Symbol
        for col in names(mirror_matrix)
            replace!(mirror_matrix[!, col], :low => :high, :high => :low)
        end
    else
        for col in names(mirror_matrix)
            replace!(mirror_matrix[!, col], 1 => -1, -1 => 1)
        end
    end
    append!(df, mirror_matrix)
end
