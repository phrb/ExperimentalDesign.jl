function build_linear_formula(factors::Array{Symbol, 1})
    linear_formula = Expr(:call)
    linear_formula.args = vcat(:+, factors)
    return Formula(:y, linear_formula)
end

function get_model_variables(formula::Formula)
    variables = Array{Any, 1}()
    push!(variables, Terms(formula).terms...)
    return variables
end

"""
    scale_orthogonal!(design::Array{Float64, 2}, factors::Array{T, 1}) where T <: Any

Orthogonally scale and center factors of a design using design and factor
limits.

# Examples

```jldoctest
julia> using ExperimentalDesign, DataStructures, StatsModels

julia> A = [5. 2. -1.; 5. 3. 0.; -5. 1. -2.]
3×3 Array{Float64,2}:
  5.0  2.0  -1.0
  5.0  3.0   0.0
 -5.0  1.0  -2.0

julia> factors = OrderedDict([(:f1, [-5., 0., 5.]), (:f2, [1., 2., 3.]), (:f3, [-2., -1., 0.])])
DataStructures.OrderedDict{Symbol,Array{Float64,1}} with 3 entries:
  :f1 => [-5.0, 0.0, 5.0]
  :f2 => [1.0, 2.0, 3.0]
  :f3 => [-2.0, -1.0, 0.0]

julia> scale_orthogonal!(A, collect(values(factors)))
3×3 Array{Float64,2}:
  1.0   0.0   0.0
  1.0   1.0   1.0
 -1.0  -1.0  -1.0

julia> A
3×3 Array{Float64,2}:
  1.0   0.0   0.0
  1.0   1.0   1.0
 -1.0  -1.0  -1.0

```

# Formula

For a design ``D_{n,p}`` with ``n`` experiments or rows and ``p`` factors
or columns, scale each factor ``\\mathbf{x}_i`` to ``\\mathbf{x}_{i}^{s}``
according to:

```math
\\mathbf{x}_{i}^{s} = \\dfrac{\\mathbf{x}_i - \\bar{M}}{\\bar{M}_{def}}
```

Where ``\\mathbf{x}_{def}`` is the factor defined in the `factors` parameter and:

```math
\\bar{M} = (max(\\mathbf{x}_i) + min(\\mathbf{x}_i)) / 2
```

```math
\\bar{M}_{def} = (max(\\mathbf{x}_{def}) - min(\\mathbf{x}_{def})) / 2
```
"""
function scale_orthogonal!(design::Array{Float64, 2},
                           factors::Array{T, 1}) where T <: Any
    for i = 1:size(design, 2)
        design_range = (max(design[:, i]...) + min(design[:, i]...)) / 2
        factor_range = (max(factors[i]...) - min(factors[i]...)) / 2
        design[:, i] = (design[:, i] .- design_range) ./ factor_range
    end

    return design
end

"""
    scale_boxdraper_encoding!(design::Array{Float64, 2},
                              factors::Array{T, 1};
                              scale_denominator = true) where T <: Any

Scale factors of a design using the Box and Draper's coding convention from
"Response Surfaces, Mixtures, and Ridge Analyses".

# Examples

With denominator scaling:

```jldoctest
julia> using ExperimentalDesign, DataStructures, StatsModels

julia> A = float(plackett_burman(4))
4×3 Array{Float64,2}:
  1.0   1.0   1.0
  1.0  -1.0  -1.0
 -1.0  -1.0   1.0
 -1.0   1.0  -1.0

julia> factors = OrderedDict([(:f1, [-1., 1.]), (:f2, [-1., 1.]), (:f3, [-1., 1.])])
DataStructures.OrderedDict{Symbol,Array{Float64,1}} with 3 entries:
  :f1 => [-1.0, 1.0]
  :f2 => [-1.0, 1.0]
  :f3 => [-1.0, 1.0]

julia> scale_boxdraper_encoding!(A, collect(values(factors)))
4×3 Array{Float64,2}:
  1.0   1.0   1.0
  1.0  -1.0  -1.0
 -1.0  -1.0   1.0
 -1.0   1.0  -1.0

julia> all([isapprox(1.0, sqrt(sum(A[:, i] .^ 2.0) / size(A, 1))) for i = 1:size(A, 2)])
true

```

Without denominator scaling:

```jldoctest
julia> using ExperimentalDesign, DataStructures, StatsModels

julia> A = float(plackett_burman(4))
4×3 Array{Float64,2}:
  1.0   1.0   1.0
  1.0  -1.0  -1.0
 -1.0  -1.0   1.0
 -1.0   1.0  -1.0

julia> factors = OrderedDict([(:f1, [-1., 1.]), (:f2, [-1., 1.]), (:f3, [-1., 1.])])
DataStructures.OrderedDict{Symbol,Array{Float64,1}} with 3 entries:
  :f1 => [-1.0, 1.0]
  :f2 => [-1.0, 1.0]
  :f3 => [-1.0, 1.0]

julia> scale_boxdraper_encoding!(A, collect(values(factors)), scale_denominator = false)
4×3 Array{Float64,2}:
  0.5   0.5   0.5
  0.5  -0.5  -0.5
 -0.5  -0.5   0.5
 -0.5   0.5  -0.5

julia> all([isapprox(1.0, sqrt(sum(A[:, i] .^ 2.0))) for i = 1:size(A, 2)])
true

```

# Formula

For a design ``D_{n,p}`` with ``n`` experiments or rows and ``p`` factors
or columns, scale each factor ``\\mathbf{x}_i`` to ``\\mathbf{x}_{i}^{s}``
according to:

```math
\\mathbf{x}_{i}^{s} = \\dfrac{\\mathbf{x}_i - \\bar{\\mathbf{x}}_i}{S_i}
```

Where ``\\bar{\\mathbf{x}}_i`` is mean of factor definition values in the
`factors` parameter and:

```math
S_{i}^{2} = \\dfrac{1}{n} \\sum\\limits_{j = 1}^{n}{(x_{ij} - \\bar{\\mathbf{x}}_i)^{2}}
```

If we pass `scale_denominator = false`, ``S_i`` becomes:

```math
S_{i}^{2} = \\sum\\limits_{j = 1}^{n}{(x_{ij} - \\bar{\\mathbf{x}}_i)^{2}}
```
"""
function scale_boxdraper_encoding!(design::Array{Float64, 2},
                                   factors::Array{T, 1};
                                   scale_denominator = true) where T <: Any
    for i = 1:size(design, 2)
        factor_mean = mean(factors[i])
        denominator = sum((design[:, i] .- factor_mean) .^ 2.0)

        if scale_denominator
            denominator /= size(design, 1)
        end

        denominator = sqrt(denominator)

        numerator = design[:, i] .- factor_mean

        if !iszero(denominator) && !iszero(numerator)
            numerator ./= denominator
        end

        design[:, i] = numerator
    end

    return design
end

"""
    get_expanded_values(factors::OrderedDict)

Return a mapping from an original level of a factor to a column
in an expanded design. Does not expand factors.

# Examples

```jldoctest
julia> using ExperimentalDesign, DataStructures, StatsModels

julia> A = OrderedDict(:A => [1., 2.], :B => [1, 2, 3, 4], :C => ["A", "B"], :D => [1.2, 2.2, 3.4], :E => [:a, :b, :c, :d])
DataStructures.OrderedDict{Symbol,Any} with 5 entries:
  :A => [1.0, 2.0]
  :B => [1, 2, 3, 4]
  :C => String["A", "B"]
  :D => [1.2, 2.2, 3.4]
  :E => Symbol[:a, :b, :c, :d]

julia> get_expanded_values(A)
DataStructures.OrderedDict{Any,Any} with 2 entries:
  :C => DataStructures.OrderedDict{Any,Symbol}("B"=>:C_1)
  :E => DataStructures.OrderedDict{Any,Symbol}(:b=>:E_1,:c=>:E_2,:d=>:E_3)

```
"""
function get_expanded_values(factors::OrderedDict)
    expanded_values = OrderedDict()
    for factor in factors
        if !(eltype(factor[2]) <: Real)
            expanded_values[factor[1]] = OrderedDict{Any, Symbol}()
            for l = 2:length(factor[2])
                expanded_values[factor[1]][factor[2][l]] = Symbol(factor[1], "_$(l - 1)")
            end
        end
    end
    return expanded_values
end

function expand_design(design::Array{T, 2}, factors::OrderedDict) where T <: Any
    small_design = DataFrame(design)
    rename!(small_design, OrderedDict(zip(names(small_design), keys(factors))))

    exp_values = get_expanded_values(factors)
    expanded_design = DataFrame()

    for name in names(small_design)
        if name in keys(exp_values)
            for key in keys(exp_values[name])
                expanded_design[exp_values[name][key]] = [small_design[i, name] == key ? 1.0 : 0. for i = 1:length(small_design[name])]
            end
        else
            expanded_design[name] = small_design[name]
        end
    end

    return expanded_design
end

"""
    expand_factors(factors::OrderedDict)

Return an `OrderedDict` with true factors expanded to `[0., 1.]` 2-level numerical factors.

For a true factor with ``n`` levels, creates ``n - 1`` 2-level factors. Only
one of these new factors can be at level `1.` at each row of a design, encoding
each level of the original factor. The first level of the original factor is
encoded by all new factors being at level `0.`.

# Examples
```jldoctest
julia> using ExperimentalDesign, DataStructures, StatsModels

julia> A = OrderedDict(:A => [1., 2.], :B => [1, 2, 3, 4], :C => ["A", "B"], :D => [1.2, 2.2, 3.4], :E => [:a, :b, :c, :d])
DataStructures.OrderedDict{Symbol,Any} with 5 entries:
  :A => [1.0, 2.0]
  :B => [1, 2, 3, 4]
  :C => String["A", "B"]
  :D => [1.2, 2.2, 3.4]
  :E => Symbol[:a, :b, :c, :d]

julia> expand_factors(A)
DataStructures.OrderedDict{Symbol,Any} with 7 entries:
  :A   => [1.0, 2.0]
  :B   => [1, 2, 3, 4]
  :C_1 => [0.0, 1.0]
  :D   => [1.2, 2.2, 3.4]
  :E_1 => [0.0, 1.0]
  :E_2 => [0.0, 1.0]
  :E_3 => [0.0, 1.0]

```
"""
function expand_factors(factors::OrderedDict)
    expanded_factors = OrderedDict{Symbol, Any}()
    for factor in factors
        if !(eltype(factor[2]) <: Real)
            for l = 1:(length(factor[2]) - 1)
                expanded_factors[Symbol(factor[1], "_$l")] = [0., 1.]
            end
        else
            expanded_factors[factor[1]] = factor[2]
        end
    end
    return expanded_factors
end

"""
    generate_model_matrix(formula::Formula,
                          design::Array{Float64, 2},
                          factors::OrderedDict;
                          scale::Function = scale_boxdraper_encoding!)

Generate a `DataFrame` with a scaled model matrix for a given formula, design and factors.

Assumes that `formula` is a linear relationship between all the factors in `factors`.

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

```
"""
function generate_model_matrix(formula::Formula,
                               input_design::DataFrame,
                               factors::OrderedDict;
                               scale::Function = scale_boxdraper_encoding!)
    design = deepcopy(input_design)
    contrasts = []

    for factor in keys(factors)
        if !(eltype(factors[factor]) <: Real)
            push!(contrasts,
                  (factor, ContrastsMatrix(DummyCoding(), factors[factor])))
        end
    end

    design[:y] = ones(size(design, 1))

    model_frame = ModelFrame(formula, design, contrasts = Dict(contrasts))

    return ModelMatrix(model_frame).m
end
