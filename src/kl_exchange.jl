using LinearAlgebra, StatsModels, ExperimentalDesign, Distributions, DataFrames

"""
$(TYPEDSIGNATURES)

```jldoctest
julia> 1 + 1
2
```
"""
function d_criterion(model_matrix;
                     tolerance = 1e-9)
    det(((model_matrix' * model_matrix) + (I * tolerance)) /
        size(model_matrix, 1)) ^
        (1 / size(model_matrix, 2))
end

"""
$(TYPEDSIGNATURES)

Builds optimum designs using the KL exchange algorithm, as described by Atkinson
*et al.*.   The ideia is  to iteratively swap  ``K`` design elements  with ``L``
candidate elements,  in order to  *maximize* an optimality  criterion.  Although
any criteria based  on information matrices could be used  here, the KL exchange
algorithm     leverages    determinant     properties     to    optimize     the
[`d_criterion`](@ref).

> Atkinson, A., Donev, A., & Tobias, R. (2007). Optimum experimental designs, with
> SAS (Vol. 34). Oxford University Press, Chapter 12.

```jldoctest
julia> candidates = FullFactorial(fill([-1, 0, 1], 5));

julia> candidates_formula = ConstantTerm(0) ~ sum(Term.(Symbol.(names(candidates.matrix))));

julia> candidates_formula = FormulaTerm(candidates_formula.lhs,
                                        candidates_formula.rhs +
                                        (@formula 0 ~ factor3 ^ 2).rhs);

julia> selected_rows = kl_exchange(candidates_formula,
                                   candidates.matrix,
                                   seed_design_size = 2,
                                   experiments = 11,
                                   design_k = 11,
                                   candidates_l = size(candidates.matrix, 1) - 11);

julia> d_criterion(selected_rows)
0.730476820204009

julia> candidates.matrix[selected_rows.indices[1], :]
11×5 DataFrames.DataFrame
│ Row │ factor1 │ factor2 │ factor3 │ factor4 │ factor5 │
│     │ Int64   │ Int64   │ Int64   │ Int64   │ Int64   │
├─────┼─────────┼─────────┼─────────┼─────────┼─────────┤
│ 1   │ -1      │ 1       │ -1      │ 1       │ 1       │
│ 2   │ 1       │ -1      │ 1       │ 1       │ 1       │
│ 3   │ 1       │ 1       │ 0       │ -1      │ 1       │
│ 4   │ -1      │ 1       │ 1       │ -1      │ 1       │
│ 5   │ -1      │ -1      │ 0       │ 1       │ 1       │
│ 6   │ 1       │ 1       │ 1       │ 1       │ -1      │
│ 7   │ -1      │ -1      │ 1       │ -1      │ -1      │
│ 8   │ 1       │ -1      │ 0       │ 1       │ -1      │
│ 9   │ -1      │ 1       │ -1      │ 1       │ -1      │
│ 10  │ -1      │ 1       │ 0       │ -1      │ -1      │
│ 11  │ 1       │ -1      │ -1      │ -1      │ -1      │
```

```jldoctest
julia> n_candidates = 3000;

julia> n_factors = 30;

julia> design_generator = DesignDistribution(Distributions.Uniform(0, 1), n_factors);

julia> candidates = rand(design_generator, n_candidates);

julia> selected_rows = kl_exchange(candidates.formula,
                                   candidates.matrix,
                                   experiments = n_factors + 2,
                                   design_k = n_factors,
                                   candidates_l = round(Int, (n_candidates - n_factors) / 2));

julia> d_criterion(selected_rows)
0.07904918814259465

julia> candidates.matrix[selected_rows.indices[1], :]
32×30 DataFrames.DataFrame. Omitted printing of 24 columns
│ Row │ factor1    │ factor2   │ factor3   │ factor4   │ factor5   │ factor6   │
│     │ Float64    │ Float64   │ Float64   │ Float64   │ Float64   │ Float64   │
├─────┼────────────┼───────────┼───────────┼───────────┼───────────┼───────────┤
│ 1   │ 0.832487   │ 0.875035  │ 0.433788  │ 0.805738  │ 0.481517  │ 0.382364  │
│ 2   │ 0.761498   │ 0.0398811 │ 0.998741  │ 0.632286  │ 0.29121   │ 0.55497   │
│ 3   │ 0.639793   │ 0.49432   │ 0.678795  │ 0.685274  │ 0.269203  │ 0.896397  │
│ 4   │ 0.00370418 │ 0.5139    │ 0.0630861 │ 0.175531  │ 0.570244  │ 0.60512   │
│ 5   │ 0.698514   │ 0.207593  │ 0.766977  │ 0.0719112 │ 0.665055  │ 0.499193  │
│ 6   │ 0.412744   │ 0.785263  │ 0.432861  │ 0.909976  │ 0.642545  │ 0.832013  │
│ 7   │ 0.353518   │ 0.116685  │ 0.140352  │ 0.929095  │ 0.0484366 │ 0.838524  │
⋮
│ 25  │ 0.577343   │ 0.790529  │ 0.763994  │ 0.548131  │ 0.402066  │ 0.293336  │
│ 26  │ 0.0941214  │ 0.52891   │ 0.293733  │ 0.92224   │ 0.982713  │ 0.929744  │
│ 27  │ 0.0256128  │ 0.723072  │ 0.23927   │ 0.993221  │ 0.521335  │ 0.407319  │
│ 28  │ 0.221631   │ 0.375189  │ 0.67251   │ 0.658613  │ 0.272448  │ 0.230501  │
│ 29  │ 0.843885   │ 0.375497  │ 0.116069  │ 0.989937  │ 0.935444  │ 0.0466424 │
│ 30  │ 0.181611   │ 0.956062  │ 0.98997   │ 0.80352   │ 0.994557  │ 0.203797  │
│ 31  │ 0.738407   │ 0.25753   │ 0.929704  │ 0.407498  │ 0.0934688 │ 0.99359   │
│ 32  │ 0.228975   │ 0.81272   │ 0.803107  │ 0.931826  │ 0.598371  │ 0.588823  │
```
"""
function kl_exchange(formula::FormulaTerm,
                     candidate_set::DataFrame;
                     tolerance::Float64 = 1e-9,
                     seed_design_size::Int = 2,
                     max_iterations::Int = 1000,
                     experiments::Int,
                     design_k::Int,
                     candidates_l::Int)
    model_matrix = ModelMatrix(ModelFrame(formula, candidate_set)).m

    seed = view(model_matrix,
                rand(1:size(model_matrix, 1), seed_design_size),
                :)

    l_variances = zeros(Int, candidates_l)
    k_variances = zeros(Int, design_k)

    design_information = zeros(size(model_matrix, 2),
                               size(model_matrix, 2))

    variances = zeros(size(model_matrix, 1) - size(seed, 1))

    # Seed Iterations
    for i = 1:(experiments - seed_design_size)
        seed_information = view(design_information, 1:size(seed, 2), 1:size(seed, 2))
        seed_information .= inv((seed' * seed) + (I * tolerance))

        candidates = view(model_matrix,
                          setdiff(1:size(model_matrix, 1), seed.indices[1]),
                          :)

        variances_candidates = view(variances, 1:size(candidates, 1))
        variances_candidates .= dot.(eachrow(candidates * seed_information),
                                     eachrow(candidates))

        l_variances .= partialsortperm(variances_candidates, 1:candidates_l, rev = true)

        seed = view(model_matrix,
                    vcat(candidates.indices[1][findmax(l_variances)[2]],
                         seed.indices[1]),
                    :)

    end

    design = seed
    variances_design = zeros(size(design, 1))
    variances_candidates = view(variances, 1:(size(model_matrix, 1) - size(design, 1)))

    kl_deltas = zeros(design_k, candidates_l)

    # Exchange Steps
    for i = 1:max_iterations
        design_information .= inv((design' * design) + (I * tolerance))
        candidates = view(model_matrix,
                          setdiff(1:size(model_matrix, 1), design.indices[1]),
                          :)

        variances_design .= dot.(eachrow(design * design_information),
                                 eachrow(design))
        k_variances .= partialsortperm(variances_design, 1:design_k, rev = false)

        variances_candidates .= dot.(eachrow(candidates * design_information),
                                     eachrow(candidates))
        l_variances .= partialsortperm(variances_candidates, 1:candidates_l, rev = true)

        k_swappable = view(model_matrix,
                           [design.indices[1][k] for k in k_variances],
                           :)

        l_swappable = view(model_matrix,
                           [candidates.indices[1][l] for l in l_variances],
                           :)

        kl_deltas .= k_swappable * design_information * l_swappable'

        swap = ((-1, -1), -1)

        for k = 1:design_k
            for l = 1:candidates_l
                delta = ((1 - variances_design[k_variances[k]]) *
                         (1 + variances_candidates[l_variances[l]])) +
                         (kl_deltas[k, l] ^ 2)

                if delta > swap[2]
                    swap = ((k, l), delta)
                end
            end
        end

        if swap[2] - 1.0 <= tolerance
            break
        end

        design = view(model_matrix,
                      vcat(deleteat!(copy(design.indices[1]),
                                     k_variances[swap[1][1]]),
                           candidates.indices[1][l_variances[swap[1][2]]]),
                      :)
    end

    return design
end

function test()
    c = 300000
    n = 50

    println("Allocating memory")
    design_generator = DesignDistribution(Distributions.Uniform(0, 1), n)
    @time candidates = rand(design_generator, c)

    println("Design created, calling exchange")
    println(design_generator.formula)

    @time b = kl_exchange(design_generator.formula,
                          candidates,
                          experiments = n + 2,
                          design_k = n - 3,
                          candidates_l = round(Int, (c - n) / 2))
    return b
end
