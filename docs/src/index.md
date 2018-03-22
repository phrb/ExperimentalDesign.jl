# ExperimentalDesigns.jl Documentation

```@docs
scale_orthogonal!(design::Array{Float64, 2}, factors::Array{T, 1}) where T <: Any

scale_boxdraper_encoding!(design::Array{Float64, 2},
                          factors::Array{T, 1};
                          scale_denominator = false) where T <: Any

generate_model_matrix(formula::DataFrames.Formula,
                      design::Array{Float64, 2},
                      factors::Array{T, 1};
                      scale::Function = scale_orthogonal!) where T <: Any

d_optimality(model_matrix::Array{Float64, 2})

d_efficiency_lower_bound(model_matrix::Array{Float64, 2})
```
