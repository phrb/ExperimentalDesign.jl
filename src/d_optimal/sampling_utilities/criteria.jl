function get_prediction_variances(model_matrix::Array{Float64, 2})
    information_matrix = model_matrix' * model_matrix

    if det(information_matrix) != 0.0
        dispersion_matrix = inv(information_matrix)
        rows              = size(dispersion_matrix, 1)

        prediction_variances = [dispersion_matrix[i, :]' * dispersion_matrix *
                                dispersion_matrix[i, :] for i = 1:rows]

        return prediction_variances
    else
        return 0.0
    end
end

"""
    d_optimality(model_matrix::Array{Float64, 2})

Compute the ``D``-optimality of a given design's model matrix.

# Examples

```jldoctest
julia> using ExperimentalDesign, DataStructures, StatsModels

julia> A = plackett_burman(4)
4×3 Array{Int64,2}:
  1   1   1
  1  -1  -1
 -1  -1   1
 -1   1  -1

julia> factors = OrderedDict([(:f1, [-1., 1.]), (:f2, [-1., 1.]), (:f3, [-1., 1.])])
DataStructures.OrderedDict{Symbol,Array{Float64,1}} with 3 entries:
  :f1 => [-1.0, 1.0]
  :f2 => [-1.0, 1.0]
  :f3 => [-1.0, 1.0]

julia> m = generate_model_matrix(@formula(y ~ f1 + f2 + f3), float(A), factors)
4×4 Array{Float64,2}:
 1.0   1.0   1.0   1.0
 1.0   1.0  -1.0  -1.0
 1.0  -1.0  -1.0   1.0
 1.0  -1.0   1.0  -1.0

julia> d_optimality(m)
256.0

```

# Formula

The ``D``-optimality of a design is the determinant of the information matrix
``\\mathbf{X}^{\\prime}\\mathbf{X}``, where ``\\mathbf{X}`` is the model matrix of
a design.

"""
function d_optimality(model_matrix::Array{Float64, 2})
    det_information_matrix = det(model_matrix' * model_matrix)

    return det_information_matrix < 0.0 ? 0.0 : det_information_matrix
end

"""
    d_efficiency_lower_bound(model_matrix::Array{Float64, 2})

Compute a lower bound for the ``D``-efficiency of a given design's model matrix
according to Castillo's "Process Optimization : A Statistical Approach".

```jldoctest
julia> using ExperimentalDesign, DataStructures, StatsModels

julia> A = plackett_burman(4)
4×3 Array{Int64,2}:
  1   1   1
  1  -1  -1
 -1  -1   1
 -1   1  -1

julia> factors = OrderedDict([(:f1, [-1., 1.]), (:f2, [-1., 1.]), (:f3, [-1., 1.])])
DataStructures.OrderedDict{Symbol,Array{Float64,1}} with 3 entries:
  :f1 => [-1.0, 1.0]
  :f2 => [-1.0, 1.0]
  :f3 => [-1.0, 1.0]

julia> m = generate_model_matrix(@formula(y ~ f1 + f2 + f3), float(A), factors)
4×4 Array{Float64,2}:
 1.0   1.0   1.0   1.0
 1.0   1.0  -1.0  -1.0
 1.0  -1.0  -1.0   1.0
 1.0  -1.0   1.0  -1.0

julia> d_efficiency_lower_bound(m)
1.0

```

# Formula

For a design ``A_{n,p}`` with ``n`` experiments or rows, ``p`` factors or
columns, and model matrix ``\\mathbf{X}``, the lower bound for the
``D``-efficiency of ``A`` is:

```math
D_{eff}^{(L)} = \\dfrac{|\\mathbf{X}^{\\prime}\\mathbf{X}|^{1/p}}{n}
```
"""
function d_efficiency_lower_bound(model_matrix::Array{Float64, 2})
    return ^(d_optimality(model_matrix), 1 / size(model_matrix, 2)) /
           size(model_matrix, 1)
end

function a_optimality(model_matrix::Array{Float64, 2})
    information_matrix = model_matrix' * model_matrix

    if det(information_matrix) != 0.0
        return trace(inv(information_matrix)) / size(model_matrix, 2)
    else
        return 0.0
    end

end

function v_optimality(model_matrix::Array{Float64, 2})
    prediction_variances = get_prediction_variances(model_matrix)
    rows                 = size(model_matrix, 1)

    return sum(prediction_variances) / rows
end

function g_optimality(model_matrix::Array{Float64, 2})
    prediction_variances = get_prediction_variances(model_matrix)
    return max(prediction_variances...)
end

function g_efficiency(model_matrix::Array{Float64, 2})
    prediction_variances = get_prediction_variances(model_matrix)
    max_variance         = max(prediction_variances...)

    g_e = size(model_matrix, 2) / max_variance

    if g_e == Inf
        return 0.0
    else
        return g_e
    end
end

function d_efficiency_lower_bound_algdesign(model_matrix::Array{Float64, 2})
    g_e   = g_efficiency(model_matrix)
    d_elb = exp(1 - (1 / g_e))

    if d_elb == Inf
        return 0.0
    else
        return d_elb
    end
end

function condition_number(model_matrix::Array{Float64, 2})
    condition_number = cond(model_matrix)

    if condition_number == Inf
        return 0.0
    else
        return condition_number
    end
end
